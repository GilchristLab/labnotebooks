library(testthat)
library(ribModel)

context("RFP Model")

#test_that("RFP Model testing simulated versus actual accuracy", {
  # Skip unless manually run or changed
  #if (F)
#    skip("RFP Model testing is optional.")
  
  #####################
  ### Initial Setup ###
  #####################

  # Test with Jeremy's simulated genome data
  fileName = file.path("UnitTestingData", "testJeremyFiles", "JeremySimulatedRFPData.csv")
  fileTable = file.path("UnitTestingData", "testJeremyFiles", "codonTranslationRates.csv")
  fileTruePhiValues = file.path("UnitTestingData", "testJeremyFiles", "RFPPhiValues.csv")
  fileTrueAlphaValues = file.path("UnitTestingData", "testJeremyFiles", "RFPAlphaValues.csv")
  fileTrueLambdaPrimeValues = file.path("UnitTestingData", "testJeremyFiles", "RFPLambdaPrimeValues.csv")
  
  # Ensure the input files exist.
  test_that("file exists: JeremySimulatedRFPData.csv", {
    expect_equal(file.exists(fileName), T)
  })
  
  test_that("file exists: codonTranslationRates.csv", {
    expect_equal(file.exists(fileTable), T)
  })
  
  test_that("file exists: RFPPhiValues.csv", {
    expect_equal(file.exists(fileTruePhiValues), T)
  })
  
  test_that("file exists: RFPAlphaValues.csv", {
    expect_equal(file.exists(fileTrueAlphaValues), T)
  })
  
  test_that("file exists: RFPLambdaPrimeValues.csv", {
    expect_equal(file.exists(fileTrueLambdaPrimeValues), T)
  })
  
  genome <- initializeGenomeObject(file = fileName, FALSE)
  
  sphi_init <- c(2)
  numMixtures <- 1
  mixDef <- "allUnique"
  geneAssignment <- c(rep(1, length(genome))) 
  parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE, mixture.definition = mixDef)
  #parameter <- initializeParameterObject(model="RFP", restart.file="30restartFile.rst")
  
  samples <- 10
  thinning <- 10
  adaptiveWidth <- 10
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                               est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
  
  model <- initializeModelObject(parameter, "RFP")
  setRestartSettings(mcmc, "restartJeremyFile.rst", adaptiveWidth, TRUE)
  
  outFile = file.path("UnitTestingOut", "testRFPJeremyLog10.txt")
  
  sink(outFile)
  system.time(
    runMCMC(mcmc, genome, model, 8)
  )
  sink()
  
  #########################################################################################
  ### Output File 1: MCMC, Loglikelihood, Mixture Probability, Mphi, Sphi, Expected Phi ###
  #########################################################################################
  
  # plots different aspects of trace
  trace <- parameter$getTraceObject()
  writeParameterObject(parameter, file = file.path("UnitTestingOut", "RFPJeremyObject.Rdat"))
  writeMCMCObject(mcmc, file = file.path("UnitTestingOut", "MCMCJeremyObject.Rdat"))
  
  
  pdf(file.path("UnitTestingOut", "RFP_Jeremy_allUnique_startCSP_startPhi_adaptSphi_true.pdf"))
  plot(mcmc, main = "MCMC Trace") #plots the whole loglikelihood trace
  
  
  # take a subset of the trace values for the logliklihood trace and plot them.
  # The primary reason for doing this is the "jump" that throws the scale of the graph
  # at the beginning is removed by taking out the beginning values.
  loglik.trace <- mcmc$getLogLikelihoodTrace()
  
  # The multiplier (currently 0.7) determines how much of the beginning trace is eliminated.
  # With 0.7, we are taking only the last 30% of the trace.
  percentTraceIgnore <- 0.7
  
  start <- length(loglik.trace) * percentTraceIgnore 
  
  logL <- logL <- mean(loglik.trace[start:length(loglik.trace)]) #get the mean for the subset
  plot(loglik.trace[start:length(loglik.trace)], type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)")
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  
  
  plot(trace, what = "MixtureProbability", main = "Mixture Probability")
  plot(trace, what = "Mphi", main = "Mphi")
  plot(trace, what = "Sphi", main = "Sphi")
  plot(trace, what = "ExpectedPhi", main = "Expected Phi")
  
  
  #loglik.trace <- mcmc$getLogLikelihoodTrace()
  #acf(loglik.trace)
  #acf(loglik.trace[start:length(loglik.trace)])
  dev.off()
  
  #################################
  ### Output File 2: CSP Traces ###
  #################################
  
  pdf(file.path("UnitTestingOut", "RFP_CSP_Values_Jeremy.pdf"), width = 11, height = 20)
  #plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
  plot(trace, what = "Alpha", mixture = 1)
  plot(trace, what = "LambdaPrime", mixture = 1)
  plot(trace, what = "MeanWaitingTime", mixture = 1)
  plot(trace, what = "VarWaitingTime", mixture = 1)
  dev.off()
  
  #####################
  ### Output File 3 ###
  #####################
  
  pdf(file.path("UnitTestingOut", "RFPJeremyConfidenceIntervalsForAlphaAndLambdaPrime.pdf"))

  #eventually this will need loop over all categories if there are multiple mixtures
  mixtureElement <- 1
  proposal <- FALSE
  codonList <- codons()
  # The true values of Alpha and Lambda Prime do not have the last 3 codons in our codon table, so exclude them.
  minCodonList <- codonList[-c(62,63,64)]
  minCodonLength <- length(minCodonList)
  
  alphaList <- numeric(minCodonLength)
  lambdaPrimeList <- numeric (minCodonLength)
  waitingTimes <- numeric(minCodonLength)
  alpha.ci <- matrix(0, ncol=2, nrow=minCodonLength)
  lambdaPrime.ci <- matrix(0, ncol=2, nrow=minCodonLength)
  psiList <- numeric(length(genome))
  ids <- numeric(length(genome))
  
  waitRates <- numeric(minCodonLength)
  for (i in 1:minCodonLength)
  {
    codon <- codonList[i]
    alphaList[i] <- parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 0, FALSE)
    alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
    alpha.ci[i,] <- quantile(alphaTrace[(samples * percentTraceIgnore):samples], probs = c(0.025,0.975))
      #parameter$getCodonSpecificQuantile(mixtureElement, samples * percentTraceIgnore, codon, 0, probs = c(0.025,0.975), FALSE)
      #quantile(alphaTrace[(samples * percentTraceIgnore):samples], probs = c(0.025,0.975))

    lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 1, FALSE)
    lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
    lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * percentTraceIgnore):samples], probs = c(0.025,0.975))
      #parameter$getCodonSpecificQuantile(mixtureElement, samples * percentTraceIgnore, codon, 1, probs = c(0.025,0.975), FALSE)
      #quantile(lambdaPrimeTrace[(samples * percentTraceIgnore):samples], probs = c(0.025,0.975))
    
    waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
    waitRates[i] <- (1.0/waitingTimes[i])
  }

  for (geneIndex in 1:length(genome))
  {
    psiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * percentTraceIgnore, geneIndex, 1)

    g <- genome$getGeneByIndex(geneIndex, F)
    ids[geneIndex] <- g$id
  }
  
  #Plot confidence intervals for alpha and lambda prime
  plot(NULL, NULL, xlim=range(1:minCodonLength, na.rm = T), ylim=range(alpha.ci),
       main = "Confidence Intervals for Alpha Parameter", xlab = "Codons",
       ylab = "Estimated values", axes=F)
  confidenceInterval.plot(x = 1:minCodonLength, y = alphaList, sd.y = alpha.ci)
  axis(2)
  axis(1, tck = 0.02, labels = codonList[1:minCodonLength], at=1:minCodonLength, las=2, cex.axis=.6)

  plot(NULL, NULL, xlim=range(1:minCodonLength, na.rm = T), ylim=range(lambdaPrime.ci),
       main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons",
       ylab = "Estimated values", axes=F)
  confidenceInterval.plot(x = 1:minCodonLength, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
  axis(2)
  axis(1, tck = 0.02, labels = codonList[1:minCodonLength], at=1:minCodonLength, las=2, cex.axis=.6)
  
  
  # correlation between RFPModel and Pop's wait rates
  # load Pop's data
  X <- read.csv(fileTable)
  X <- X[order(X[,1]) , ]
  XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)


  Y <- data.frame(codonList[-c(62,63,64)], waitRates)
  colnames(Y) <- c("Codon", "PausingTimeRates")
  Y <- Y[order(Y[,1]) , ]


  plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]),
       main = "Correlation Between Pop and RFP Model Pausing Time Rates", xlab = "Pop's Rates", ylab = "RFP's Rates")
  upper.panel.plot(XM[,2], Y[,2])

  # correlation between RFPModel WAIT RATES (inverse) and Pop's wait rates
  Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
  colnames(Y) <- c("Codon", "WaitingTimeRates")
  Y <- Y[order(Y[,1]) , ]


  plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]),
       main = "Correlation Between Pop and RFP Model Waiting Times", xlab = "Pop's Rates", ylab = "RFP's Waiting Times")
  upper.panel.plot(XM[,2], Y[,2])

  dev.off()
  
  #########################
  ### Estimated vs True ###
  #########################
  
  pdf(file.path("UnitTestingOut", "RFPJeremyEstimatedVsTrue.pdf"))
  
  # This is a very roughly written bit of code to more-or-less "only" print the points
  # as well as two arbitrary-looking lines for convenience.
  # TODO: Make this better/cleaner.
  
  ###########
  ### Phi ###
  ###########
  X <- read.csv(fileTruePhiValues)
  X <- X[order(X[,2]) , ]
  #XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  #fX <- XM[,2]
  
  # Vector of doubles
  fX <- X[,2]
  
  estimValues <- numeric(length(genome))
  for (geneIndex in 1:length(genome))
  {
    estimValues[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * percentTraceIgnore, geneIndex, mixtureElement)
    #tail(trace$getSynthesisRateTraceForGene(geneIndex), n=1)
  }
  
  Y <- data.frame(estimValues)
  colnames(Y) <- c("Jeremy's Phi Values")
  # Vector of doubles
  fY <- Y[order(Y[,1]) , ]
  
  plot(NULL, NULL, log = "xy", xlim=range(fX, na.rm = T), ylim=range(fY), 
       main = "Jeremy's Estimated Phi Values vs True Values", xlab = "True Values", ylab = "Jeremy's Phi Values")
  upper.panel.plot(fX, fY)
  
  #############
  ### Alpha ###
  #############
  X <- read.csv(fileTrueAlphaValues)
  X <- X[order(X[,2]) , ]
  #XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  #fX <- XM[,2]
  
  # Vector of doubles
  fX <- X[,2]
  
  estimValues <- numeric(length(genome))
  for (codonIndex in 1:length(minCodonList))
  {
    estimValues[codonIndex] <- parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 0, FALSE)
      #alphaList[codonIndex]
      #parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 0, FALSE)
  }
  
  Y <- data.frame(estimValues)
  colnames(Y) <- c("Jeremy's Alpha Values")
  # Vector of doubles
  fY <- Y[order(Y[,1]) , ]
  
  plot(NULL, NULL, xlim=range(fX, na.rm = T), ylim=range(fY), 
       main = "Jeremy's Estimated Alpha Values vs True Values", xlab = "True Values", ylab = "Jeremy's Alpha Values")
  upper.panel.plot(fX, fY)
  
  ####################
  ### Lambda Prime ###
  ####################
  X <- read.csv(fileTrueLambdaPrimeValues)
  X <- X[order(X[,2]) , ]
  #XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  #fX <- XM[,2]
  
  # Vector of doubles
  fX <- X[,2]
  
  estimValues <- numeric(length(genome))
  for (codonIndex in 1:length(minCodonList))
  {
    estimValues[codonIndex] <- parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 1, FALSE)
      #lambdaPrimeList[codonIndex]
      #parameter$getCodonSpecificPosteriorMean(mixtureElement, samples * percentTraceIgnore, codon, 1, FALSE)
  }
  
  Y <- data.frame(estimValues)
  colnames(Y) <- c("Jeremy's Lambda Prime Values")
  # Vector of doubles
  fY <- Y[order(Y[,1]) , ]
  
  plot(NULL, NULL, xlim=range(fX, na.rm = T), ylim=range(fY), 
       main = "Jeremy's Estimated Lambda Prime Values vs True Values", xlab = "True Values", ylab = "Jeremy's Lambda Prime Values")
  upper.panel.plot(fX, fY)
  
  dev.off()
  
  #################
  ### CSV Files ###
  #################
  # Will write csv files based off posterior for alpha, lambda prime, and psi
  m <- matrix(c(minCodonList, alphaList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "Alpha")
  write.csv(m, file.path("UnitTestingOut", "RFPJeremyAlphaValues.csv"), quote = F, row.names = F)
  
  
  m <- matrix(c(minCodonList, lambdaPrimeList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "LambdaPrime")
  write.table(m, file.path("UnitTestingOut", "RFPJeremyLambdaPrimeValues.csv"), quote = F, row.names = F)
  
  
  m <- matrix(c(ids, psiList, psiList), ncol = 3, byrow = FALSE)
  colnames(m) <- c("Gene", "PsiValue", "PsiValue")
  write.table(m, file.path("UnitTestingOut", "RFPJeremyPsiValues.csv"), quote = F, row.names = F)
#})
