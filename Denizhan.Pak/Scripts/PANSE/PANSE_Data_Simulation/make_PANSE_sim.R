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
  rfp_file = "orderedPopPADataRand500.csv"
  phi_file = "orderedRandGeneIDPhiMean.csv"
  
  # Test with Jeremy's simulated genome data
  fileName = file.path("../Input/PopData/", rfp_file)
  fileTable = file.path("../Input/PopData/", phi_file)
  fileTruePhiValues = file.path("../Input/PopData/", phi_file)
  fileTrueAlphaValues = file.path("../Input/PopData/", "JeremyRFPAlphaValues.csv")
  fileTrueLambdaPrimeValues = file.path("../Input/PopData/", "JeremyRFPLambdaPrimeValues.csv")
  fileTrueNSEValues = file.path("../Input/Stochastic_Simulation/", "simNSEMay_4.csv")
  
  
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
  
  test_that("file exists: simNSEMay_4.csv", {
    expect_equal(file.exists(fileTrueNSEValues), T)
  })
  
  genome <- initializeGenomeObject(file = fileName, fasta = FALSE)
  
  phiValues <- read.table(file = fileTable, sep = ",", header = TRUE)
  phiMean <- phiValues[,2]
  sphi_init <- c(1)
  numMixtures <- 1
  mixDef <- "allUnique"
  
  geneAssignment <- c(rep(1, length(genome))) #
  parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, 
                                         initial.expression.values = phiMean, model= "PANSE", split.serine = TRUE, mixture.definition = mixDef)
  parameter$initMutationSelectionCategories(c(fileTrueAlphaValues), 1, "Alpha")
  parameter$initMutationSelectionCategories(c(fileTrueLambdaPrimeValues), 1, "LambdaPrime")
  parameter$initMutationSelectionCategories(c(fileTrueNSEValues), 1, "NSERate")
  
  model <- initializeModelObject(parameter, "PANSE", rfp.count.column = 1)
  model$simulateGenome(genome)
  genome$writeRFPData("Pop_Sim_With_Errors.csv", T)
  
