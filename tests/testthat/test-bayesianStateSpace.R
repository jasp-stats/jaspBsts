
# NOTE: the tests assume an English locale, otherwise things with dates fail
oldLocale <- Sys.getlocale("LC_TIME")
if (oldLocale != "en_US.UTF-8") {
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  on.exit(Sys.setlocale("LC_TIME", oldLocale))
}

os <- Sys.info()[['sysname']]


dateTimeVarNames <- c("dateYear", "dateWeek", "dateDay", "timeHour", "timeMinute", "timeSecond")
#### tests on windows ####


test_that("Model Summary table results match (AR - manual)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$autoregressiveComponent <- TRUE
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.411857605442053, 1.05018027614687, 0.542967391044393, 0.811701682421661
                                   ))
  }

})

test_that("Model Summary table results match (AR - automatic)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$autoregressiveComponent <- TRUE
  options$lagSelectionMethod <- "auto"
  options$maxLags <- 4
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.535488123460269, 1.05773803700731, 0.539404673068745, 0.72136258779103
                                   ))
  }

})


test_that("Model Summary table results match (local level)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$localLevelComponent <- TRUE
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0432055353493266, 1.09380783381569, 0.520147149047757, 1.08103597020621
                                   ))
  }

})


test_that("Model Summary table results match (local linear trend)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$localLinearTrendComponent <- TRUE
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0777667265276332, 1.20104889746718, 0.423311739669594, 1.09879731436253
                                   ))
  }

})

test_that("Model Summary table results match (semi-local linear trend)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$localLevelComponent <- TRUE
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0432055353493266, 1.09380783381569, 0.520147149047757, 1.08103597020621
                                   ))
  }

})

test_that("Model Summary table results match (seasonal)", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$localLevelComponent <- TRUE
  options$seasonalities <- list(list(normalPriorMean = 0,
                                     numbers = 2,
                                     name = "season",
                                     inverseGammaPriorN = 0.01,
                                     duration = 1,
                                     normalPriorSd="",
                                     inverseGammaPriorSd=""
  ))
  options$samples <- 10

  set.seed(1)
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0394752039178554, 1.12753310072425, 0.490267461631722, 1.07910143584074
                                   ))
  }

})

options <- jaspTools::analysisOptions("bayesianStateSpace")
options$dependent <- "contNormal"
options$covariates <- c("contcor1", "contcor2")
options$time <- "dateYear"
options$posteriorSummaryTable <- TRUE
options$localLevelComponent <- TRUE
options$samples <- 10
options$expectedPredictors <- 3
options$modelTerms <- list(list(components = "contcor1", isNuisance = FALSE), list(
  components = "contcor2", isNuisance = FALSE))
set.seed(1)


test_that("Model Summary table results match (local level - covariates)", {
  skip_on_os(c("mac","linux"))

  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.0446411563379834, 1.07342874539745, 0.538811093349419, 1.03451899091434
                                   ))
  }

})

test_that("Posterior Summary of Coefficients table results match (local level - covariates)", {
  skip_on_os(c("mac","linux"))
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list("<unicode><unicode><unicode>", "contcor1", 0.0710805371191776,
                                        0.314637544172676, 1, 1, "", 0.375148071900377, "<unicode><unicode><unicode>",
                                        "contcor2", -0.287511723102937, -0.11088682615349, 1, 1, "",
                                        0.0421351582659255, 0, "(Intercept)", 0, 0, 0, 1, "", 0))

  }

})


options <- jaspTools::analysisOptions("bayesianStateSpace")
options$dependent <- "contNormal"
options$fixedFactors <- c("facGender", "facExperim")
options$time <- "dateYear"
options$posteriorSummaryTable <- TRUE
options$localLevelComponent <- TRUE
options$samples <- 10
options$expectedPredictors <- 3
options$modelTerms <- list(list(components = "facGender", isNuisance = FALSE), list(
  components = "facExperim", isNuisance = FALSE))
set.seed(1)


test_that("Model Summary table results match (local level - fixedFactors)", {
  skip_on_os(c("mac","linux"))
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.0579395417976514, 1.06599469439761, 0.546337972019031, 1.02729362023781
                                   ))
  }

})

test_that("Posterior Summary of Coefficients table results match (local level - fixedFactors)", {
  skip_on_os(c("mac","linux"))
  for (i in 1:6){
    options$time <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list("<unicode><unicode><unicode>", "facGenderm", 0.253231187396135,
                                        0.609442717859108, 1, 1, "", 0.638874711991056, "<unicode><unicode><unicode>",
                                        "facExperimexperimental", -0.222080709900353, 0.0399657905691665,
                                        1, 1, "", 0.271697756233028, 0, "(Intercept)", 0, 0, 0, 1, "",
                                        0))

  }

})

options$aggregatedStatesPlot <- TRUE
options$componentStatesPlot <- TRUE
options$forecastErrorPlot<- TRUE
options$residualPlot <- TRUE
options$controlChartPlot <- TRUE
options$probalisticControlPlot <- TRUE
for (i in 1:6){
  options$time <- dateTimeVarNames[i]

  results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

  test_that("State plots matches", {
    skip_on_os(c("mac","linux"))
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsComponentStatePlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("component-states-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsAggregatedStatePlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("aggregated-state-",i))
  })
  test_that("Residual Plots match", {
    skip_on_os(c("mac","linux"))
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsResidualPlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("residual-plot-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsForecastErrorPlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("forecast-error-plot-",i))
  })
  test_that("Control chart plots match", {
    skip_on_os(c("mac","linux"))
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotThreshold"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("threshold-plot-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotProbability"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("probability-plot-",i))
  })
}

test_that("Prediction plot matches", {
  skip_on_os(c("mac","linux"))
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$samples <- 10
  options$localLevelComponent <- TRUE
  options$time <- "dateYear"
  options$predictionHorizon <- 12


  results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
  plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsPredictionPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "prediction-plot")
})








#### tests on mac/linux ####

  test_that("Model Summary table results match (AR - manual)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$autoregressiveComponent <- TRUE
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(0.314607509367699, 1.05068512793848, 0.543107992387219, 0.876243631114363
                                     ))
    }

  })

  test_that("Model Summary table results match (AR - automatic)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$autoregressiveComponent <- TRUE
    options$lagSelectionMethod <- "auto"
    options$maxLags <- 4
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(0.549871515593311, 1.05769993725513, 0.539387282554917, 0.710106440387526
                                     ))
    }

  })


  test_that("Model Summary table results match (local level)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$localLevelComponent <- TRUE
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(-0.0762491443534843, 1.09400650414813, 0.519902709350384, 1.09802344438273
                                     ))
    }

  })


  test_that("Model Summary table results match (local linear trend)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$localLinearTrendComponent <- TRUE
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(-0.110766726762034, 1.20050404586317, 0.423811661472609, 1.11549244792232
                                     ))
    }

  })

  test_that("Model Summary table results match (semi-local linear trend)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$localLevelComponent <- TRUE
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(-0.0762491443534843, 1.09400650414813, 0.519902709350384, 1.09802344438273
                                     ))
    }

  })

  test_that("Model Summary table results match (seasonal)", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$localLevelComponent <- TRUE
    options$seasonalities <- list(list(normalPriorMean = 0,
                                       numbers = 2,
                                       name = "season",
                                       inverseGammaPriorN = 0.01,
                                       duration = 1,
                                       normalPriorSd="",
                                       inverseGammaPriorSd=""
    ))
    options$samples <- 10

    set.seed(1)
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(-0.0781263744430147, 1.12728326564074, 0.490397326469068, 1.09898063195421
                                     ))
    }

  })

  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$covariates <- c("contcor1", "contcor2")
  options$time <- "dateYear"
  options$posteriorSummaryTable <- TRUE
  options$localLevelComponent <- TRUE
  options$samples <- 10
  options$expectedPredictors <- 3
  options$modelTerms <- list(list(components = "contcor1", isNuisance = FALSE), list(
    components = "contcor2", isNuisance = FALSE))
  set.seed(1)


  test_that("Model Summary table results match (local level - covariates)", {
    skip_on_os("windows")

    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(0.00213853172231659, 1.0755654786762, 0.53451984213854, 1.05728076187172
                                     ))
    }

  })

  test_that("Posterior Summary of Coefficients table results match (local level - covariates)", {
    skip_on_os("windows")
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list("<unicode>", "contcor1", 0.0437460818434989, 0.275790341778563,
                                          1, 1, 0.142009518230305, 0.434946375191193, "<unicode>", "contcor2",
                                          -0.321932103241681, -0.173635046880455, 1, 1, 0.105996027757027,
                                          -0.0310099257617435, 0, "(Intercept)", 0, 0, 0, 1, 0, 0))

    }

  })


  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$fixedFactors <- c("facGender", "facExperim")
  options$time <- "dateYear"
  options$posteriorSummaryTable <- TRUE
  options$localLevelComponent <- TRUE
  options$samples <- 10
  options$expectedPredictors <- 3
  options$modelTerms <- list(list(components = "facGender", isNuisance = FALSE), list(
    components = "facExperim", isNuisance = FALSE))
  set.seed(1)


  test_that("Model Summary table results match (local level - fixedFactors)", {
    skip_on_os("windows")

    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list(0.0153869242271224, 1.06857850992189, 0.542811453062196, 1.05023866481121
                                     ))
    }

  })

  test_that("Posterior Summary of Coefficients table results match (local level - fixedFactors)", {
    skip_on_os("windows")
    for (i in 1:6){
      options$time <- dateTimeVarNames[i]
      results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
      table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
      jaspTools::expect_equal_tables(table,
                                     list("<unicode>", "facGenderm", 0.136625688072959, 0.461928614855287,
                                          1, 1, 0.220891318762768, 0.742348181051067, "<unicode>", "facExperimexperimental",
                                          -0.275546362761004, -0.0501265437197341, 1, 1, 0.163978897116324,
                                          0.172523964528349, 0, "(Intercept)", 0, 0, 0, 1, 0, 0))

    }

  })

  options$aggregatedStatesPlot <- TRUE
  options$componentStatesPlot <- TRUE
  options$forecastErrorPlot<- TRUE
  options$residualPlot <- TRUE
  options$controlChartPlot <- TRUE
  options$probalisticControlPlot <- TRUE
  # the 3rd value in dateTimeVarNames causes errors when unit test performed on Github (day)
  # possibly due to different region settings on operating systems, removed for now
  for (i in c(1:2,4:6)){
    options$time <- dateTimeVarNames[i]

    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

    test_that("State plots matches", {
      skip_on_os("windows")
      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsComponentStatePlot"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("component-states-unix-",i))

      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsAggregatedStatePlot"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("aggregated-state-unix-",i))
    })
    test_that("Residual Plots match", {
      skip_on_os("windows")
      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsResidualPlot"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("residual-plot-unix-",i))

      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsForecastErrorPlot"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("forecast-error-plot-unix-",i))
    })
    test_that("Control chart plots match", {
      skip_on_os("windows")
      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotThreshold"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("threshold-plot-unix-",i))

      plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotProbability"]][["data"]]
      testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
      jaspTools::expect_equal_plots(testPlot, paste0("probability-plot-unix-",i))
    })
  }

  test_that("Prediction plot matches", {
    skip_on_os("windows")
    options <- jaspTools::analysisOptions("bayesianStateSpace")
    options$dependent <- "contNormal"
    options$samples <- 10
    options$localLevelComponent <- TRUE
    options$time <- "dateYear"
    options$predictionHorizon <- 12


    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsPredictionPlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, "prediction-plot-unix-")
  })






