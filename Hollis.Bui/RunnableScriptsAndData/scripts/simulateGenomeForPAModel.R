library(ribModel)
rm(list=ls())

genome <- initializeGenomeObject(file.path("..", "JeremyData", "s288c.genome.fasta"))
alpha.file <- file.path("..", "JeremyData", "JeremyRFPAlphaValues.csv")
lambdaPrime.file <- file.path("..", "JeremyData", "JeremyRFPLambdaPrimeValues.csv")
phi.file <- file.path("..", "JeremyData", "JeremyRFPPhiValues.csv")
genome.out.file <- file.path("..", "HollisSimulatedGenomes", "HollisSimulatedGenomeNew.fasta")

sphi_init <- c(2)
numMixtures <- 1
mixDef <- "allUnique"

geneAssignment <- rep(0, length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model= "PA", split.serine = TRUE,
                                       mixture.definition = mixDef)
model <- initializeModelObject(parameter, "PA")

parameter$initMutationSelectionCategories(c(alpha.file), 1, "Alpha")
parameter$initMutationSelectionCategories(c(lambdaPrime.file), 1, "LambdaPrime")

phi.vals <- parameter$readPhiValues(phi.file)
parameter$initializeSynthesisRateByList(phi.vals)

model$simulateGenome(genome)

genome$writeFasta(genome.out.file, TRUE)
