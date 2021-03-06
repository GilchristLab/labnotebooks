library(ribModel)
library(VGAM)
genome <- ribModel::initializeGenomeObject("./RibModelDev/data/rfp/rfp.counts.by.codon.and.gene.GSE63789.wt.csv", fasta = FALSE)
parameter <- ribModel::initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)), model = "RFP")
mcmc <- ribModel::initializeMCMCObject(samples = 500, thinning = 10, adaptive.width=50)
model <- ribModel::initializeModelObject(parameter = parameter, model = "RFP")
ribModel::runMCMC(mcmc = mcmc, genome = genome, model = model)
