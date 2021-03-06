library(testthat)
library(AnaCoDa)

#test_that("PANSE Model testing simulated versus actual accuracy", {
  # Skip unless manually run or changed
  #if (F)
#    skip("PANSE Model testing is optional.")
  
  #####################
  ### Initial Setup ###
  #####################
  set.seed(500)
  rfp_file = "simulated_rfp_file_500_genes.csv"
  phi_file = "simulated_phi_file_500_genes.csv"
  
  # Test with Jeremy's simulated genome data
  fileName = file.path("../Input/Stochastic_Simulation/", rfp_file)
  fileTable = file.path("../Input/Stochastic_Simulation/", phi_file)
  fileTruePhiValues = file.path("../Input/Stochastic_Simulation/", phi_file)
  fileTrueAlphaValues = file.path("../Input/Stochastic_Simulation/", "simAlphaJanuary.csv")
  fileTrueLambdaPrimeValues = file.path("../Input/Stochastic_Simulation/", "simLambdaPrimeJanuary.csv")
  
  # Ensure the input files exist.
  test_that("file exists: simulated_rfp_file_750_genes.csv", {
    expect_equal(file.exists(fileName), T)
  })
  
  test_that("file exists: simulated_phi_file_750_genes.csv", {
    expect_equal(file.exists(fileTable), T)
  })
  
  test_that("file exists: simulated_phi_file_750_genes.csv", {
    expect_equal(file.exists(fileTruePhiValues), T)
  })
  
  test_that("file exists: RFPAlphaValues.csv", {
    expect_equal(file.exists(fileTrueAlphaValues), T)
  })
  
  test_that("file exists: RFPLambdaPrimeValues.csv", {
    expect_equal(file.exists(fileTrueLambdaPrimeValues), T)
  })
  
  
  genome <- initializeGenomeObject(file = fileName, fasta = FALSE)
  
  phiValues <- read.table(file = fileTable, sep = ",", header = TRUE)
  phiMean <- phiValues[,2]
  sphi_init <- c(1)
  numMixtures <- 1
  mixDef <- "allUnique"
  
  geneAssignment <- c(rep(1, length(genome))) 
  parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, 
                                         initial.expression.values = phiMean, model= "PA", split.serine = TRUE, mixture.definition = mixDef)
  #parameter <- initializeParameterObject(model="RFP", restart.file="30restartFile.rst")
  parameter$initMutationSelectionCategories(c(fileTrueAlphaValues), 1, "Alpha")
  parameter$initMutationSelectionCategories(c(fileTrueLambdaPrimeValues), 1, "LambdaPrime")
  model <- initializeModelObject(parameter, "PA", rfp.count.column = 1)
  
  samples <- 500
  thinning <- 20
  adaptiveWidth <- 10
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning, adaptive.width = adaptiveWidth, 
                               est.expression=F, est.csp=T, est.hyper=F)
  
  logFileString <- paste0("runPAModel2018Log", samples)
  outFile = file.path("./Out", logFileString)
  setRestartSettings(mcmc, "./Out/Restart/restartPAAugust2018File.rst", adaptiveWidth, TRUE)
  
  
  sink(outFile)
  system.time(
    runMCMC(mcmc, genome, model)
  )
  sink()
  
  
  #########################################################################################
  ### Output File 1: MCMC, Loglikelihood, Mixture Probability, Mphi, Sphi, Expected Phi ###
  #########################################################################################
  
  # plots different aspects of trace
  trace <- parameter$getTraceObject()
  
  writeParameterObject(parameter, file = file.path("Out", "RFPJeremyObject.Rdat"))
  writeMCMCObject(mcmc, file = file.path("Out", "MCMCJeremyObject.Rdat"))
  
  
  pdf(file.path("./Out", "RFP_Jeremy_allUnique_startCSP_startPhi_adaptSphi_true.pdf"))
  plot(mcmc, main = "MCMC Trace") #plots the whole loglikelihood trace
  
  
   # take a subset of the trace values for the logliklihood trace and plot them.
   # The primary reason for doing this is the "jump" that throws the scale of the graph
   # at the beginning is removed by taking out the beginning values.
   loglik.trace <- mcmc$getLogLikelihoodTrace()
   start <- length(loglik.trace) * 0.7
   # the multiplier (currently 0.7) determines how much of the beginning trace is eliminated.
  
   logL <- mean(loglik.trace[start:length(loglik.trace)]) #get the mean for the subset
   plot(loglik.trace[start:length(loglik.trace)], type="l", main=paste("logL:", logL), xlab="Sample", ylab="log(Likelihood)")
   grid (NULL,NULL, lty = 6, col = "cornsilk2")
  
   acceptanceRate.trace <- mcmc$getLogLikelihoodTrace()
   accRate <- mean(acceptanceRate.trace) #get the mean for the subset
   plot(acceptanceRate.trace, type="l", main=paste("Acceptance Rate:", accRate), xlab="Sample", ylab="log(Likelihood)")
   grid (NULL,NULL, lty = 6, col = "cornsilk2")
  
   plot(trace, what = "MixtureProbability", main = "Mixture Probability")
   plot(trace, what = "Mphi", main = "Mphi")
   plot(trace, what = "Sphi", main = "Sphi")
   plot(trace, what = "ExpectedPhi", main = "Expected Phi")
   loglik.trace <- mcmc$getLogLikelihoodTrace()
   acf(loglik.trace)
   acf(loglik.trace[start:length(loglik.trace)])
   dev.off()
  
  #################################
  ### Output File 2: CSP Traces ###
  #################################
  
  pdf(file.path("Out", "RFP_CSP_Values_Jeremy.pdf"), width = 11, height = 20)
  #plot(trace, what = "Expression", geneIndex = 905) #used to make sure gene stabalized, not really needed now
  plot(trace, what = "Alpha", mixture = 1)
  plot(trace, what = "LambdaPrime", mixture = 1)
  plot(trace, what = "MeanWaitingTime", mixture = 1)
  plot(trace, what = "VarWaitingTime", mixture = 1)
  plot(trace, what = "AcceptanceRatio", mixture = 1)
  
  dev.off()
  
  #####################
  ### Output File 3 ###
  #####################
  
  pdf(file.path("Out", "RFPJeremyConfidenceIntervalsForAlphaAndLambdaPrime.pdf"))

  #eventually this will need loop over all categories if there are multiple mixtures
  cat <- 1
  proposal <- FALSE
  alphaList <- numeric (61)
  lambdaPrimeList <- numeric (61)
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
    waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
    waitRates[i] <- (1.0/waitingTimes[i])
  }

  for (geneIndex in 1:length(genome))
  {
    psiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples * 0.5, geneIndex, 1)

    g <- genome$getGeneByIndex(geneIndex, F)
    ids[geneIndex] <- g$id
  }
  
  #Plot confidence intervals for alpha and lambda prime
  plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(alpha.ci),
       main = "Confidence Intervals for Alpha Parameter", xlab = "Codons",
       ylab = "Estimated values", axes=F)
  confidenceIntervalplot(x = 1:61, y = alphaList, sd.y = alpha.ci)
  axis(2)
  axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)

  plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(lambdaPrime.ci),
       main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons",
       ylab = "Estimated values", axes=F)
  confidenceIntervalplot(x = 1:61, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
  axis(2)
  axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
  
  
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
  upperpanelplot(XM[,2], Y[,2])

  # correlation between RFPModel WAIT RATES (inverse) and Pop's wait rates
  Y <- data.frame(codonList[-c(62,63,64)], waitingTimes)
  colnames(Y) <- c("Codon", "WaitingTimeRates")
  Y <- Y[order(Y[,1]) , ]


  plot(NULL, NULL, xlim=range(XM[,2], na.rm = T), ylim=range(Y[,2]),
       main = "Correlation Between Pop and RFP Model Waiting Times", xlab = "Pop's Rates", ylab = "RFP's Waiting Times")
  upperpanelplot(XM[,2], Y[,2])

  #########################
  ### Estimated vs True ###
  #########################
  
  # This is a very roughly written bit of code to more-or-less "only" print the points
  # as well as two arbitrary-looking lines for convenience.
  # TODO: Make this better/cleaner.
  
  # Phi
  X <- read.csv(fileTruePhiValues)
  X <- X[order(X[,2]) , ]
  XM <- matrix(c(X[,1], X[,2]), ncol = 2, byrow = FALSE)
  # Vector of doubles
  fX <- XM[,2]
  
  estimValues <- numeric(length(genome))
  for (geneIndex in 1:length(genome))
  {
    estimValues[geneIndex] <- tail(trace$getSynthesisRateTraceForGene(geneIndex), n=1)
  }
  
  Y <- data.frame(estimValues)
  colnames(Y) <- c("Jeremy's RFP Rates")
  # Vector of doubles
  fY <- Y[order(Y[,1]) , ]
  
  plot(NULL, NULL, log = "xy", xlim=range(fX, na.rm = T), ylim=range(fY), 
       main = "Jeremy's Estimated Phi Values vs True Values", xlab = "True Values", ylab = "Jeremy's RFP Rates")
  upper.panel.plot(fX, fY)
  
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
  
  dev.off()
  
  #################
  ### CSV Files ###
  #################
  # Will write csv files based off posterior for alpha, lambda prime, and psi
  m <- matrix(c(codonList[-c(62,63,64)], alphaList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "Alpha")
  write.csv(m, file.path("TestingOut", "RFPJeremyAlphaValues.csv"), quote = F, row.names = F)
  
  
  m <- matrix(c(codonList[-c(62,63,64)], lambdaPrimeList), ncol = 2, byrow = FALSE)
  colnames(m) <- c("Codon", "LambdaPrime")
  write.table(m, file.path("TestingOut", "RFPJeremyLambdaPrimeValues.csv"), quote = F, row.names = F)
  
  
  m <- matrix(c(ids, psiList, psiList), ncol = 3, byrow = FALSE)
  colnames(m) <- c("Gene", "PsiValue", "PsiValue")
  write.table(m, file.path("TestingOut", "RFPJeremyPsiValues.csv"), quote = F, row.names = F)
#})

  
