library(ribModel)
rm(list=ls())

genome <- initializeGenomeObject(file.path("HollisTestingData", "s288c.genome.fasta"))
alpha.file <- file.path("HollisTestingData", "RFPAlphaValues.csv")
lambdaPrime.file <- file.path("HollisTestingData", "RFPLambdaPrimeValues.csv")
phi.file <- file.path("HollisTestingData", "RFPPhiValues.csv")
genome.out.file <- file.path("HollisTestingOut", "HollisSimulatedGenome.fasta")

sphi_init <- c(2)
numMixtures <- 1
mixDef <- "allUnique"

geneAssignment <- rep(0, length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "RFP", split.serine = TRUE,
                                       mixture.definition = mixDef)
model <- initializeModelObject(parameter, "RFP")

parameter$initMutationSelectionCategories(c(alpha.file), 1, "Alpha")
parameter$initMutationSelectionCategories(c(lambdaPrime.file), 1, "LambdaPrime")

phi.vals <- parameter$readPhiValues(phi.file)
parameter$initializeSynthesisRateByList(phi.vals)

model$simulateGenome(genome)

genome$writeFasta(genome.out.file, TRUE)
