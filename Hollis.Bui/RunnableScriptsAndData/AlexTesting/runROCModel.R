library(ribModel) #change here
rm(list=ls())
#read genome

with.phi <- FALSE

if (with.phi) {
  genome <- initializeGenomeObject(file = "RibModelDev/data/twoMixtures/simulatedAllUniqueR.fasta", expression.file = "RibModelDev/data/twoMixtures/simulatedAllUniqueR_phi.csv")
} else {
  genome <- initializeGenomeObject(file = "Ecoli_K12_MG1655_sigpep_nt.fasta")
}
cat("Genome loaded\n")
#initialize parameter object
sphi_init <- 1
numMixtures <- 1
mixDef <- "allUnique"
size <- length(genome)
index <- c(1:size)
geneAssignment <- rep(1, length(genome))

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)

#parameter <- initializeParameterObject(restart.file = "2330restartFile.rst")
# initialize MCMC object
samples <-100
thinning <- 50
adaptiveWidth <- 50
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                     est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)

# get model object
#parameter <- loadParameterObject("signal_peptide_genes.Rdat")
model <- initializeModelObject(parameter, "ROC", with.phi)

setRestartSettings(mcmc, "restartFile.rst", adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 4)
)

#plots log likelihood trace, possibly other mcmc diagnostics in the future
plot(mcmc)

getCSPEstimates(parameter,"selection_ecoli_sigpep.csv","Selection",1,samples*0.1)
getCSPEstimates(parameter,"mutation_ecoli_sigpep.csv","Mutation",1,samples*0.1)


mixtureAssignment <- getMixtureAssignmentEstimate(parameter,size,samples*0.1)
expressionValues <- getExpressionEstimatesForMixture(parameter,size,mixtureAssignment,samples*0.1)

expressionValues <- log10(expressionValues)
write(expressionValues,file="ecoli_gene_expression_just_sigpep.txt",ncolumns = 1)
#writeParameterObject(parameter,"signal_peptide_genes.Rdat")
# plots different aspects of trace
trace <- parameter$getTraceObject()
pdf("ecoli_k12_mg1655_just_sigpep.pdf")
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")
dev.off()
plot(trace, what = "Expression", geneIndex = 63)
pdf("trace_values_ecoli_just_sigpep.pdf", width = 11, height = 12)
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
plot(model, genome, samples = samples*0.1, mixture = 1, main = "E.coli Codon Usage Plot")
dev.off()
