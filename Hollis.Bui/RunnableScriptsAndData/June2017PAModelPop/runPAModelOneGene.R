library(testthat)
library(ribModel)

  #####################
  ### Initial Setup ###
  #####################

  # Define input and output directories here.
  testingIn = "TestingIn"
  testingOut = "OutputOneGene"

  # Define: Pop Data subset, true Pop codon translation rates, Gilchrist Phi means
  fileName = file.path(testingIn, "PopPADataOneGene.csv")
  fileTable = file.path(testingIn, "codonTranslationRates.csv")
  filePhiValues = file.path(testingIn, "preprocessedFiles", "geneIDPhiMean.csv")
  
  #fileTrueAlphaValues = file.path(testingIn, "JeremyRFPAlphaValues.csv")
  #fileTrueLambdaPrimeValues = file.path(testingIn, "JeremyRFPLambdaPrimeValues.csv")
  
  # Ensure the input files exist.
  test_that("file exists: PopPADataOneGene.csv", {
    expect_equal(file.exists(fileName), T)
  })
  
  test_that("file exists: codonTranslationRates.csv", {
    expect_equal(file.exists(fileTable), T)
  })
  
  test_that("file exists: geneIDPhiMean.csv", {
    expect_equal(file.exists(filePhiValues), T)
  })
  
  # test_that("file exists: JeremyRFPAlphaValues.csv", {
  #   expect_equal(file.exists(fileTrueAlphaValues), T)
  # })
  # 
  # test_that("file exists: JeremyRFPLambdaPrimeValues.csv", {
  #   expect_equal(file.exists(fileTrueLambdaPrimeValues), T)
  # })
  
  genome <- initializeGenomeObject(file = fileName, fasta = FALSE)
  
  phiValues <- read.csv(file = filePhiValues, header = TRUE)
  phiMean <- phiValues[,2]
  sphi_init <- c(2)
  numMixtures <- 1
  mixDef <- "allUnique"
  geneAssignment <- c(rep(1, length(genome)))
  
  parameter <- initializeParameterObject(genome=genome, sphi=sphi_init, num.mixtures=numMixtures, gene.assignment=geneAssignment, model="PA", split.serine=TRUE, mixture.definition=mixDef)
  #parameter <- initializeParameterObject(model="PA", restart.file="30restartFile.rst")

  samples <- 1000
  thinning <- 10
  adaptiveWidth <- 10 # Gilchrist comment -- seems large
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                               est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
  
  model <- initializeModelObject(parameter=parameter, model="PA", rfp.count.column=1)
  
  restartFileString = file.path(testingOut, "restartFiles", "restartPAModelScript2017File")
  setRestartSettings(mcmc, restartFileString, adaptiveWidth, TRUE)
  
  logFileString <- paste0("runPAModelScript2017Log", samples, ".txt")
  outFile = file.path(testingOut, logFileString)
  
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
  writeParameterObject(parameter, file = file.path(testingOut, "PAModelObject2017.Rdat"))
  writeMCMCObject(mcmc, file = file.path(testingOut, "MCMCPAModel2017.Rdat"))
  
  
  pdf(file.path(testingOut, "PA_Model2017_allUnique_startCSP_startPhi_adaptSphi_true.pdf"))
  plot(mcmc, main = "MCMC Trace") #plots the whole loglikelihood trace
  
  
  # take a subset of the trace values for the logliklihood trace and plot them.
  # The primary reason for doing this is the "jump" that throws the scale of the graph
  # at the beginning is removed by taking out the beginning values.
  loglik.trace <- mcmc$getLogLikelihoodTrace()
  start <- length(loglik.trace) * 1.0
  # the multiplier (currently 0.7) determines how much of the beginning trace is eliminated.

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
  
  pdf(file.path(testingOut, "CSP_Values_PA_Model2017.pdf"), width = 11, height = 20)
  #plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
  plot(trace, what = "Alpha", mixture = 1)
  plot(trace, what = "LambdaPrime", mixture = 1)
  plot(trace, what = "MeanWaitingTime", mixture = 1)
  plot(trace, what = "VarWaitingTime", mixture = 1)
  dev.off()
  
  #####################
  ### Output File 3 ###
  #####################
  
  pdf(file.path(testingOut, "PAModel2017ConfidenceIntervalsForAlphaAndLambdaPrime.pdf"))

  #eventually this will need loop over all categories if there are multiple mixtures
  cat <- 1
  proposal <- FALSE
  alphaList <- numeric(61)
  lambdaPrimeList <- numeric(61)
  waitingTimes <- numeric(61)
  alpha.ci <- matrix(0, ncol=2, nrow=61)
  lambdaPrime.ci <- matrix(0, ncol=2, nrow=61)
  psiList <- numeric(length(genome))
  ids <- numeric(length(genome))
  codonList <- codons()

  waitRates <- numeric(61)
  for (i in 1:61)
  {
    codon <- codonList[i]
    alphaList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 0, FALSE)
    alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
    alpha.ci[i,] <- quantile(alphaTrace[(samples * 0.5):samples], probs = c(0.025,0.975))

    lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 1, FALSE)
    lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
    lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
    
    waitingTimes[i] <- alphaList[i] / lambdaPrimeList[i]
    waitRates[i] <- (1.0/waitingTimes[i])
  }

  for (geneIndex in 1:length(genome))
  {
    psiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)

    g <- genome$getGeneByIndex(geneIndex, F)
    ids[geneIndex] <- g$id
  }
  
  #Plot confidence intervals for alpha and lambda prime
  # plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(alpha.ci),
  #      main = "Confidence Intervals for Alpha Parameter", xlab = "Codons",
  #      ylab = "Estimated values", axes=F)
  # confidenceInterval.plot(x = 1:61, y = alphaList, sd.y = alpha.ci)
  # axis(2)
  # axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
  # 
  # plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(lambdaPrime.ci),
  #      main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons",
  #      ylab = "Estimated values", axes=F)
  # confidenceInterval.plot(x = 1:61, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
  # axis(2)
  # axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
  
  
  # # correlation between RFPModel and Pop's wait rates
  # # load Pop's data
  # X <- read.csv(fileTable)
  # X <- X[order(X[,1]) , ]
  # XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  # 
  # 
  # Y <- data.frame(codonList[-c(62,63,64)], waitRates)
  # colnames(Y) <- c("Codon", "PausingTimeRates")
  # Y <- Y[order(Y[,1]) , ]
  # 
  # 
  # plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]),
  #      main = "Correlation Between Pop and Our PA Model Pausing Time Rates", xlab = "Pop's Rates", ylab = "PA Model's Rates")
  # upper.panel.plot(XM[,2], Y[,2])
  # 
  # # correlation between PA Model WAIT RATES (inverse) and Pop's wait rates
  # Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
  # colnames(Y) <- c("Codon", "WaitingTimeRates")
  # Y <- Y[order(Y[,1]) , ]
  # 
  # 
  # plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]),
  #      main = "Correlation Between Pop and PA Model Waiting Times", xlab = "Pop's Rates", ylab = "PA Model's Waiting Times")
  # upper.panel.plot(XM[,2], Y[,2])

  #########################
  ### Estimated vs True ###
  #########################
  
  # This is a very roughly written bit of code to more-or-less "only" print the points
  # as well as two arbitrary-looking lines for convenience.
  # TODO: Make this better/cleaner.
  
  # Phi
  #X <- read.csv(filePhiValues)
  # X <- read.table(file = filePhiValues, sep = '\t', header = TRUE)
  # X <- X[order(X[,2]) , ]
  # XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  # # Vector of doubles
  # fX <- XM[,2]
  # 
  # estimValues <- numeric(length(genome))
  # for (geneIndex in 1:length(genome))
  # {
  #   estimValues[geneIndex] <- tail(trace$getSynthesisRateTraceForGene(geneIndex), n=1)
  # }
  # 
  # Y <- data.frame(estimValues)
  # colnames(Y) <- c("PA Model Rates")
  # # Vector of doubles
  # fY <- Y[order(Y[,1]) , ]
  # 
  # plot(NULL, NULL, log = "xy", xlim=range(fX, na.rm = T), ylim=range(fY), 
  #      main = "PA Model Estimated Phi Values vs True Values", xlab = "True Values", ylab = "PA Model RFP Rates")
  # upper.panel.plot(x = fX, y = fY)
  
  # # Alpha
  # X <- read.csv(fileTrueAlphaValues)
  # X <- X[order(X[,2]) , ]
  # XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  # # Vector of doubles
  # fX <- XM[,2]
  # 
  # estimValues <- numeric(length(genome))
  # for (geneIndex in 1:length(genome))
  # {
  #   estimValues[geneIndex] <- tail(trace$getCodonSpecificParameterTraceByMixtureElement(1, codonList[geneIndex], 0, 0), n=1)
  # }
  # 
  # Y <- data.frame(estimValues)
  # colnames(Y) <- c("Jeremy's RFP Rates")
  # # Vector of doubles
  # fY <- Y[order(Y[,1]) , ]
  # 
  # plot(NULL, NULL, log = "xy", xlim=range(fX, na.rm = T), ylim=range(fY), 
  #      main = "Jeremy's Estimated Alpha Values vs True Values", xlab = "True Values", ylab = "Jeremy's RFP Rates")
  # upper.panel.plot(fX, fY)
  # ###
  
  ### TEST
  
  # plot(parameter, what = "Mutation", samples = samples)
  
  # dev.off()
  
  #################
  ### CSV Files ###
  #################
  # Will write csv files based off posterior for alpha, lambda prime, and psi
  m <- matrix(c(codonList[-c(62,63,64)], alphaList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "Alpha")
  write.csv(m, file.path(testingOut, "PAModel2017AlphaValues.csv"), quote = F, row.names = F)
  
  
  m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "LambdaPrime")
  write.csv(m, file.path(testingOut, "PAModel2017LambdaPrimeValues.csv"), quote = F, row.names = F)
  
  
  # TODO: Make a third column here for Phi values, probably.
  m <- matrix(c(ids, psiList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Gene", "PsiValue")
  write.csv(m, file.path(testingOut, "PAModel2017PsiValues.csv"), quote = F, row.names = F)
#})
