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

# This is a generated file. Don't change it

bayesianStateSpace <- function(
          data = NULL,
          version = "0.17",
          formula = NULL,
          isNuisance = NULL,
          aggregatedStatesPlot = FALSE,
          aggregatedStatesPlotCiLevel = 0.95,
          aggregatedStatesPlotObservationsShown = FALSE,
          automaticBurninProportion = 0.1,
          autoregressiveComponent = FALSE,
          burninMethod = "auto",
          componentStatesPlot = FALSE,
          controlChartPlot = FALSE,
          controlPeriod = 100,
          controlSigma = 2,
          covariates = list(),
          dependent = "",
          dynamicRegregressionComponent = FALSE,
          dynamicRegregressionLags = 0,
          expectedPredictors = 1,
          fixedFactors = list(),
          forecastErrorPlot = FALSE,
          lagSelectionMethod = "manual",
          lags = 1,
          localLevelComponent = FALSE,
          localLinearTrendComponent = FALSE,
          manualBurninAmount = 0,
          maxLags = 1,
          modelTerms = list(),
          plotHeight = 320,
          plotWidth = 480,
          posteriorSummaryCiLevel = 0.95,
          posteriorSummaryTable = FALSE,
          predictionHorizon = 0,
          probalisticControlPlot = FALSE,
          residualPlot = FALSE,
          samples = 2000,
          seasonalities = list(),
          seed = 1,
          semiLocalLinearTrendComponent = FALSE,
          showCoefMeanInc = FALSE,
          time = "",
          timeout = 120) {

   defaultArgCalls <- formals(jaspBsts::bayesianStateSpace)
   defaultArgs <- lapply(defaultArgCalls, eval)
   options <- as.list(match.call())[-1L]
   options <- lapply(options, eval)
   defaults <- setdiff(names(defaultArgs), names(options))
   options[defaults] <- defaultArgs[defaults]
   options[["data"]] <- NULL
   options[["version"]] <- NULL

   if (!is.null(formula)) {
      if (!inherits(formula, "formula")) {
         formula <- as.formula(formula)
      }
      options$formula <- jaspBase::jaspFormula(formula, data)
   }

   optionsWithFormula <- c("isNuisance", "covariates", "dependent", "fixedFactors", "modelTerms", "seasonalities", "time")
   for (name in optionsWithFormula) {
      if ((name %in% optionsWithFormula) && inherits(options[[name]], "formula")) options[[name]] = jaspBase::jaspFormula(options[[name]], data)   }

   return(jaspBase::runWrappedAnalysis("jaspBsts::bayesianStateSpace", data, options, version))
}