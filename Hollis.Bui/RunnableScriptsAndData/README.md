# Contents:

1. `scripts`: R files to run.
  1. `simulateRFP.R` -- A failing script to simulate genomes based on `JeremyData`.
  Kept for debugging purposes, as it should work but crashes in Rstudio when `model$simulateGenome`   is run.
  2. `testJeremyRestart.R` -- Script that requires restart files (not included) and the simulated genome file used to generate these restart files (included) to plot an MCMC run.
  3. `toCSV.pl` -- Perl script to simply transform Lareau data to .csv format for PANSE reading. Input file not included for size constraints.

2. `HollisSimulatedGenomes`: Simulated genomes produced by `JeremyData` and running `simulateGenome` in C++ only.
  1. `HollisSimulatedGenome` -- A genome simulated using the `JeremyData` file `s288c.genome.fasta`, `RFPAlphaValues.csv`, `RFPLambdaPrimeValues.csv`, and `RFPPhiValues.csv`, based on Dr. Gilchrist's suggestion of genome file to use.
  2. `HollisSimulatedGenome2` -- A genome simulated using the `JeremyData` file `rfp.counts.by.codon.and.gene.GSE63789.wt.csv`, and the .csv files described above. This simulation is based on what Gabriel had done in the past, where he did simulate a genome using this genome file.

3. `JeremyData` -- data files for input extracted from existing, true RFP work done by Jeremy. Or Gabriel.
  1. `JeremySimulatedRFPData.csv` -- A genome simulated by Jeremy via unknown means.
  2. `RFPAlphaValues` -- The true values of the genome described by GSE63789.
  3. `RFPLambdaPrimeValues.csv`
  4. `RFPPhiValues.csv`
  5. `codonTranslationRates.csv` -- The codon translation rates (inverse of our pausing times) described by GSE63789.
  6. `rfp.counts.by.codon.and.gene.GSE63789.wt.csv` -- The true rfp values of the genome described by GSE63789
  7. `s288c.genome.fasta` -- A fasta file describing the genome observed in GSE63789.

4. `runLogs`: Various finished runs that are archived.
  1. `7-27-16-JeremyCorrelation` -- Output from running an MCMC script with 20,000 samples, 10 thinning, 10 adaptive width, and using the input files from `JeremyData` with Jeremy's simulated genome.
  2. `testRFPModelRunsWithSamples` -- Running a certain number of sampled runs with a non-simulated genome (`rfp.counts.by.codon.and.gene.GSE3789.wt.csv`, I believe) to reveal that there's some non-useful data being outputted, and some debugging may need to be done (hence the dated runs with simulated genomes and comparing them to non-simulated data).

