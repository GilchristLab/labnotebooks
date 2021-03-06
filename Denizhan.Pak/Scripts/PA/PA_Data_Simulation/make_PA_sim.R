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
  rfp_file = "simulated_rfp_file_1000_genes.csv"
  phi_file = "simulated_phi_file_1000_genes.csv"
  
  # Test with Jeremy's simulated genome data
  fileName = file.path("../Input/Stochastic_Simulation", rfp_file)
  fileTable = file.path("../Input/Stochastic_Simulation", phi_file)
  fileTruePhiValues = file.path("../Input/Stochastic_Simulation", phi_file)
  fileTrueAlphaValues = file.path("../Input/Stochastic_Simulation", "simAlphaJanuary.csv")
  fileTrueLambdaPrimeValues = file.path("../Input/Stochastic_Simulation", "simLambdaPrimeJanuary.csv")
  
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
  
  
  genome <- initializeGenomeObject(file = fileName, fasta = FALSE, simulated = F)
  
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
  model$simulateGenome(genome)
  genome$writeRFPData("RFP_test_Jan_30.csv", T)
  
