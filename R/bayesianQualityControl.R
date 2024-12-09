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

  saveRDS(dataset, file = "C:/JASP/dataset.RDS")
  saveRDS(options, file = "C:/JASP/options.RDS")

  return()

  if(.bqcReady(options))
    .bqcFitModel(jaspResults, dataset, options)
  options <- readRDS(file = "C:/JASP/options.RDS")
  dataset <- readRDS(file = "C:/JASP/dataset.RDS")

}

.bqcReady        <- function(options) {

  return(options[["measurement"]] != "")
}
.bqcDependencies <- c(
  'advancedMcmcAdaptation', 'advancedMcmcBurnin', 'advancedMcmcChains', 'advancedMcmcSamples', 'advancedMcmcThin', 'advancedStateAggregation',
  'controlLimitsLower', 'controlLimitsUpper', 'measurement', 'time'
)

.bqcFitModel  <- function(jaspResults, dataset, options) {

  if (is.null(jaspResults[["fit"]]$object))
    return()

  # create the output container
  fitContainer <- createJaspState()
  fitContainer$dependOn(.bqcReady)
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

    # Priors on innovation scales
    # sigma_mu ~ exponential(1) -> sigma_mu ~ dexp(1)
    sigma_mu ~ dexp(1)
    sigma_log_sigma ~ dexp(1)

    # State evolution
    for (t in 2:nStates) {
      # mu[t] ~ normal(mu[t-1], sigma_mu)
      # JAGS: mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))
      mu[t] ~ dnorm(mu[t-1], 1/(sigma_mu^2))

      # log_sigma[t] ~ normal(log_sigma[t-1], sigma_log_sigma)
      # JAGS: log_sigma[t] ~ dnorm(log_sigma[t-1], 1/(sigma_log_sigma^2))
      log_sigma[t] ~ dnorm(log_sigma[t-1], 1/(sigma_log_sigma^2))
    }

    # Observation model
    for (i in 1:nObservations) {
      y[i] ~ dnorm(mu[state[i]], 1/(sigma[t]^2))
    }

    # Derived quantities (analogous to generated quantities in Stan)
    for (t in 1:nStates) {
      sigma[t] <- exp(log_sigma[t])
      val1[t] <- (USL - mu[t])/(3 * sigma[t])
      val2[t] <- (mu[t] - LSL)/(3 * sigma[t])
      # cpk[t] = min(val1[t], val2[t]) -> use ifelse in JAGS
      cpk[t] <- ifelse(val1[t] < val2[t], val1[t], val2[t])
    }
  }
  "

  # Run the JAGS model using runjags
  fit <- runjags::run.jags(
    model = jagsModel,        # The JAGS model code stored as a string
    data = jagsData,          # The data list
    monitor = c("mu", "sigma", "cpk"),  # Variables to monitor
    n.chains = options[["advancedMcmcChains"]],           # Number of MCMC chains
    adapt  = options[["advancedMcmcAdaptation"]],         # Number of adaptation steps
    burnin = options[["advancedMcmcBurnin"]],             # Number of burn-in steps
    sample = options[["advancedMcmcSamples"]]             # Number of posterior samples
  )

  # Extract and save the samples
  samples <- coda::as.mcmc(fit)
  fitContainer$object <- samples

  return()
}
.bqcPlotState <- function(jaspResults, dataset, options, parameter) {

  if (is.null(jaspResults[["fit"]]$object))
    return()

  # Create the output container
  plotContainer <- createJaspPlot(title = switch(
    parameter,
    "mu"    = "Mean",
    "sigma" = "Standard deviation",
    "cpk"   = "Cpk"))
  plotContainer$dependOn(c(.bqcReady, switch(
    parameter,
    "mu"    = "statePlotsMean",
    "sigma" = "statePlotsStandardDeviation",
    "cpk"   = "statePlotsCpk")
  ))
  plotContainer$position <- switch(
    parameter,
    "mu"    = 1,
    "sigma" = 2,
    "cpk"   = 3
  )
  jaspResults[[paste0("statePlot-")]] <- plotContainer

  # Extract the samples
  samples <- jaspResults[["fit"]]$object

  # Select the parameter
  samples <- samples[, grep(parameter, colnames(samples))]
  samplesSummary <- data.frame(
    time = 1:ncol(samples),
    mean = apply(samples, 2, mean),
    lowerCI = apply(samples, 2, quantile, 0.025),
    upperCI = apply(samples, 2, quantile, 0.975)
  )

  # Create the plot
  plot <- ggplot2::ggplot(samplesSummary, ggplot2::aes(x=time, y=mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lowerCI, ymax=upperCI), alpha=0.3, fill="skyblue") +
    ggplot2::geom_line(color="blue") +
    ggplot2::geom_point(color="blue") +
    ggplot2::ggtitle("Posterior of mu over time") +
    ggplot2::ylab(switch(
      parameter,
      "mu"    = "Mean",
      "sigma" = "Standard deviation",
      "cpk"   = "Cpk")) +
    ggplot2::xlab("State") +
    jaspGraphs::geom_rangeframe(sides = "bl") +
    jaspGraphs::themeJaspRaw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  plotContainer$object <- plot

  return()
}
