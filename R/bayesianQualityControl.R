#
# Copyright (C) 2013-2022 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Main function ----

bayesianQualityControl         <- function() return(NULL)
bayesianQualityControlInternal <- function(jaspResults, dataset, options) {

  # saveRDS(dataset, file = "C:/JASP/dataset.RDS")
  # saveRDS(options, file = "C:/JASP/options.RDS")

  if (.bqcReady(options))
    .bqcFitModel(jaspResults, dataset, options)

  if (options[["statePlotsMean"]])
    .bqcPlotState(jaspResults, dataset, options, "mu")
  if (options[["statePlotsStandardDeviation"]] && options[["scaleModel"]] == "Heterogeneous")
    .bqcPlotState(jaspResults, dataset, options, "sigma")
  if (options[["statePlotsCpk"]])
    .bqcPlotState(jaspResults, dataset, options, "cpk")

  if (options[["posteriorPlotsMean"]])
    .bqcPlotPosterior(jaspResults, dataset, options, "mu")
  if (options[["posteriorPlotsStandardDeviation"]])
    .bqcPlotPosterior(jaspResults, dataset, options, "sigma")
  if (options[["posteriorPlotsCpk"]])
    .bqcPlotPosterior(jaspResults, dataset, options, "cpk")

  if (options[["CpkThresholdPlot"]])
    .bqcPlotCpkThreshold(jaspResults, dataset, options)

  return()

  options <- readRDS(file = "C:/JASP/options.RDS")
  dataset <- readRDS(file = "C:/JASP/dataset.RDS")

}

.bqcReady        <- function(options) {

  return(options[["measurement"]] != "")
}
.bqcDependencies <- c(
  'advancedMcmcAdaptation', 'advancedMcmcBurnin', 'advancedMcmcChains', 'advancedMcmcSamples', 'advancedMcmcThin', 'advancedStateAggregation',
  'controlLimitsLower', 'controlLimitsUpper', 'measurement', 'time', 'scaleModel'
)

.bqcFitModel         <- function(jaspResults, dataset, options) {

  if (!is.null(jaspResults[["fit"]]$object))
    return()

  # create the output container
  fitContainer <- createJaspState()
  fitContainer$dependOn(.bqcDependencies)
  jaspResults[["fit"]] <- fitContainer

  ### transform the data
  # create bins each containing advancedStateAggregation observations
  dataset$binnedTime <- rep(1:ceiling(nrow(dataset)/options[["advancedStateAggregation"]]), each=options[["advancedStateAggregation"]])[1:nrow(dataset)]

  jagsData <- list(
    y             = dataset[[options[["measurement"]]]],
    nObservations = nrow(dataset),
    state   = dataset$binnedTime,
    nStates = max(dataset$binnedTime),
    LSL = options[["controlLimitsLower"]],
    USL = options[["controlLimitsUpper"]]
  )
  jagsModel <- .bqcGetJagsCode(options)

  # run with 10 progress bar ticks
  jaspBase::startProgressbar(10, gettext("Estimating Bayesian State-Space Model"))

  options[["advancedMcmcSamples"]]

  # Run the JAGS model using runjags
  fit <- runjags::run.jags(
    model = jagsModel,        # The JAGS model code stored as a string
    data = jagsData,          # The data list
    monitor = c("mu", "sigma", "cpk"),  # Variables to monitor
    n.chains = options[["advancedMcmcChains"]],           # Number of MCMC chains
    adapt  = options[["advancedMcmcAdaptation"]],         # Number of adaptation steps
    burnin = options[["advancedMcmcBurnin"]],             # Number of burn-in steps
    sample = 0,
    summarise = FALSE
  )

  for(i in 1:10) {
    jaspBase::progressbarTick()
    fit <- runjags::extend.jags(fit, adapt  = 0, sample = ceiling(options[["advancedMcmcSamples"]]/10), summarise = FALSE)
  }

  # Extract and save the samples
  samples <- coda::as.mcmc(fit)
  fitContainer$object <- samples

  return()
}
.bqcPlotState        <- function(jaspResults, dataset, options, parameter) {

  if (is.null(jaspResults[["statePlotContainer"]])) {
    statePlotContainer <- createJaspContainer(title = gettext("State Plots"))
    statePlotContainer$dependOn(.bqcDependencies)
    statePlotContainer$position <- 1
    jaspResults[["statePlotContainer"]] <- statePlotContainer
  } else{
    statePlotContainer <- jaspResults[["statePlotContainer"]]
  }

  if (is.null(statePlotContainer[[parameter]])) {
    # Create the output container
    plotContainer <- createJaspPlot(title = switch(
      parameter,
      "mu"    = gettext("Mean"),
      "sigma" = gettext("Standard deviation"),
      "cpk"   = gettext("Cpk")), width = 500, height = 400)
    plotContainer$dependOn(switch(
      parameter,
      "mu"    = "statePlotsMean",
      "sigma" = "statePlotsStandardDeviation",
      "cpk"   = "statePlotsCpk")
    )
    plotContainer$position <- switch(
      parameter,
      "mu"    = 1,
      "sigma" = 2,
      "cpk"   = 3
    )
    statePlotContainer[[parameter]] <- plotContainer

  } else {
    plotContainer <- statePlotContainer[[parameter]]
  }

  if (is.null(jaspResults[["fit"]]$object))
    return()

  # Extract the samples
  samples <- jaspResults[["fit"]]$object

  # Select the parameter
  samples <- samples[, grep(parameter, colnames(samples))]
  samplesSummary <- data.frame(
    time = 1:ncol(samples),
    mean    = apply(samples, 2, median),
    lowerCI = apply(samples, 2, quantile, 0.025),
    upperCI = apply(samples, 2, quantile, 0.975)
  )

  # Create the plot
  plot <- ggplot2::ggplot(data = samplesSummary) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(x = time, ymin = lowerCI, ymax = upperCI), alpha = 0.3) +
    jaspGraphs::geom_line(mapping = ggplot2::aes(x = time, y = mean)) +
    jaspGraphs::scale_x_continuous(name = gettext("State"), breaks = jaspGraphs::getPrettyAxisBreaks(range(samplesSummary$time))) +
    jaspGraphs::scale_y_continuous(name = switch(
      parameter,
      "mu"    = gettext("Mean"),
      "sigma" = gettext("Standard deviation"),
      "cpk"   = expression("C"["pk"])), breaks = jaspGraphs::getPrettyAxisBreaks(c(min(samplesSummary$lowerCI), max(samplesSummary$upperCI)))) +
    jaspGraphs::geom_rangeframe(sides = "bl") +
    jaspGraphs::themeJaspRaw()

  plotContainer$plotObject <- plot

  return()
}
.bqcPlotPosterior    <- function(jaspResults, dataset, options, parameter) {

  if (is.null(jaspResults[["posteriorPlotContainer"]])) {
    posteriorPlotContainer <- createJaspContainer(title = gettext("Posterior Distributions"))
    posteriorPlotContainer$dependOn(.bqcDependencies)
    posteriorPlotContainer$position <- 2
    jaspResults[["posteriorPlotContainer"]] <- posteriorPlotContainer
  } else {
    posteriorPlotContainer <- jaspResults[["posteriorPlotContainer"]]
  }

  if (is.null(posteriorPlotContainer[[parameter]])) {
    # Create the output container
    plotContainer <- createJaspPlot(title = switch(
        parameter,
        "mu"    = gettext("Mean"),
        "sigma" = gettext("Standard deviation"),
        "cpk"   = "Cpk"
      ), width = 450, height = 400)
    plotContainer$position <- switch(
      parameter,
      "mu"    = 1,
      "sigma" = 2,
      "cpk"   = 3
    )
    plotContainer$dependOn(switch(
      parameter,
      "mu"    = "posteriorPlotsMean",
      "sigma" = "posteriorPlotsStandardDeviation",
      "cpk"   = "posteriorPlotsCpk")
    )
    posteriorPlotContainer[[parameter]] <- plotContainer

  } else {
    plotContainer <- posteriorPlotContainer[[parameter]]
  }

  if (is.null(jaspResults[["fit"]]$object))
    return()

  # Extract the samples
  samples <- jaspResults[["fit"]]$object

  if (parameter == "sigma" && options[["scaleModel"]] != "Heterogeneous")
    parName <- "sigma"
  else
    parName <- sprintf("%1$s[%2$i]", parameter, options[["posteriorPlotsAtState"]])

  # Check whether state exists
  if (!parName %in% colnames(samples)) {
    plotContainer$setError(gettextf("No samples for %1$s[%2$i]", parameter, options[["posteriorPlotsAtState"]]))
    return()
  }

  # Select the parameter
  samples <- samples[,parName]

  # TODO: change in future: keep only 99.5% of central samples
  samples <- samples[order(samples)]
  samples <- samples[ceiling(0.005 * length(samples) / 2):floor(0.995 * length(samples) / 2)]

  # Create a data frame for plotting
  plotData <- data.frame(value = as.vector(samples))

  # Create the plot (histogram + density line)
  plot <- ggplot2::ggplot(data = plotData, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, fill = "gray80", color = "black") +
    ggplot2::geom_density(size = 1) +
    jaspGraphs::scale_x_continuous(
      name = switch(
        parameter,
        "mu"    = gettext("Mean"),
        "sigma" = gettext("Standard deviation"),
        "cpk"   = expression("C"["pk"]),
        parameter
      ),
      breaks = jaspGraphs::getPrettyAxisBreaks(range(plotData$value))
    ) +
    jaspGraphs::scale_y_continuous(
      name = gettext("Density"),
      breaks = jaspGraphs::getPrettyAxisBreaks
    ) +
    jaspGraphs::geom_rangeframe(sides = "bl") +
    jaspGraphs::themeJaspRaw()

  plotContainer$plotObject <- plot

  return()
}
.bqcPlotCpkThreshold <- function(jaspResults, dataset, options) {

  if (!is.null(jaspResults[["cpkThresholdPlot"]]))
    return()

  plotContainer <- createJaspPlot(title = gettext("Cpk Threshold Plot"), width = 500, height = 400)
  plotContainer$dependOn(c(.bqcDependencies, "CpkThresholdPlot", "CpkThresholdPlotThreshold"))
  plotContainer$position <- 3
  jaspResults[["cpkThresholdPlot"]] <- plotContainer

  if (is.null(jaspResults[["fit"]]$object))
    return()

  # Extract the samples
  samples <- jaspResults[["fit"]]$object

  # Select the parameter
  samples <- samples[, grep("cpk", colnames(samples))]
  samples <- samples > options[["CpkThresholdPlotThreshold"]]
  samplesSummary <- data.frame(
    time    = 1:ncol(samples),
    mean    = apply(samples, 2, mean)
  )

  # Create the plot
  plot <- ggplot2::ggplot(data = samplesSummary) +
    jaspGraphs::geom_line(mapping = ggplot2::aes(x = time, y = mean)) +
    jaspGraphs::scale_x_continuous(name = gettext("State"), breaks = jaspGraphs::getPrettyAxisBreaks(range(samplesSummary$time))) +
    jaspGraphs::scale_y_continuous(name = bquote(P("C"["pk"])>.(options[["CpkThresholdPlotThreshold"]])), breaks = jaspGraphs::getPrettyAxisBreaks(c(0, 1))) +
    jaspGraphs::geom_rangeframe(sides = "bl") +
    jaspGraphs::themeJaspRaw()

  plotContainer$plotObject <- plot

  return()
}
.bqcGetJagsCode   <- function(options) {

  if (options[["scaleModel"]] == "Heterogeneous") {
    jagsModel <- "
    model {
      # Priors for initial states
      # mu[1] ~ normal((LSL + USL)/2, (USL - LSL)/(2*1.96))
      # In JAGS: dnorm(mean, tau), tau = 1/sd^2
      # sd = (USL - LSL)/(2*1.96)
      # tau = 1 / sd^2 = (2 * 1.96 / (USL - LSL))^2
      mu[1] ~ dnorm((LSL + USL)/2, pow((2 * 1.96 / (USL - LSL)), 2))

      # log_sigma[1] ~ normal(0,1) in Stan
      # In JAGS: dnorm(0,1) means mean=0, tau=1, so sd=1
      log_sigma[1] ~ dnorm(0, 1)
      sigma[1]    <- exp(log_sigma[1])

      # Priors on innovation scales
      # sigma_mu ~ exponential(1) -> sigma_mu ~ dexp(1)
      sigma_mu        ~ dexp(1)
      sigma_log_sigma ~ dexp(1)

      # State evolution
      for (t in 2:nStates) {
        # mu[t] ~ normal(mu[t-1], sigma_mu)
        # JAGS: mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))
        mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))

        # log_sigma[t] ~ normal(log_sigma[t-1], sigma_log_sigma)
        # JAGS: log_sigma[t] ~ dnorm(log_sigma[t-1], 1/(sigma_log_sigma^2))
        log_sigma[t] ~ dnorm(log_sigma[t-1], 1/(sigma_log_sigma^2))
        sigma[t]    <- exp(log_sigma[t])
      }

      # Observation model
      for (i in 1:nObservations) {
        y[i] ~ dnorm(mu[state[i]], 1/(sigma[state[i]]^2))
      }

      # Derived quantities (analogous to generated quantities in Stan)
      for (t in 1:nStates) {
        val1[t] <- (USL - mu[t])/(3 * sigma[t])
        val2[t] <- (mu[t] - LSL)/(3 * sigma[t])
        # cpk[t] = min(val1[t], val2[t]) -> use ifelse in JAGS
        cpk[t] <- ifelse(val1[t] < val2[t], val1[t], val2[t])
      }
    }
    "
  } else if (options[["scaleModel"]] == "Homogeneous") {
    jagsModel <- "
    model {
      # Priors for initial states
      # mu[1] ~ normal((LSL + USL)/2, (USL - LSL)/(2*1.96))
      # In JAGS: dnorm(mean, tau), tau = 1/sd^2
      # sd = (USL - LSL)/(2*1.96)
      # tau = 1 / sd^2 = (2 * 1.96 / (USL - LSL))^2
      mu[1] ~ dnorm((LSL + USL)/2, pow((2 * 1.96 / (USL - LSL)), 2))

      # log_sigma[1] ~ normal(0,1) in Stan
      # In JAGS: dnorm(0,1) means mean=0, tau=1, so sd=1
      log_sigma   ~ dnorm(0, 1)
      sigma      <- exp(log_sigma)

      # Priors on innovation scales
      # sigma_mu ~ exponential(1) -> sigma_mu ~ dexp(1)
      sigma_mu        ~ dexp(1)

      # State evolution
      for (t in 2:nStates) {
        # mu[t] ~ normal(mu[t-1], sigma_mu)
        # JAGS: mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))
        mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))
      }

      # Observation model
      for (i in 1:nObservations) {
        y[i] ~ dnorm(mu[state[i]], 1/(sigma^2))
      }

      # Derived quantities (analogous to generated quantities in Stan)
      for (t in 1:nStates) {
        val1[t] <- (USL - mu[t])/(3 * sigma)
        val2[t] <- (mu[t] - LSL)/(3 * sigma)
        # cpk[t] = min(val1[t], val2[t]) -> use ifelse in JAGS
        cpk[t] <- ifelse(val1[t] < val2[t], val1[t], val2[t])
      }
    }
    "
  }

  return(jagsModel)
}
