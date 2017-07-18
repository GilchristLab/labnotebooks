# testing Input files

Contents:

0. `preprocessedFiles`: Files used to create final input files after using Perl scripts.
    1. `PopNames.out`: `PopPAData.csv` processed through `printPopNames.pl`: Just the geneIDs of the file.
    2. `PopPAData.csv`: Processed version of the Pop data with RFPcounts (aka rfp.count...) from Dr. Gilchrist's MATLAB code, with slight reordering done so it's in PA format.
    3. `PopPADataRand500`: Using `randNames.out` and `PopPAData.csv`, use `hashRandNamesPopData.pl` to create a subset of `PopPAData.csv` based on this set of random geneIDs.
    4. `geneIDPhiMean.csv`: Using data from `phiValues009670-6.tsv` and `PopNames.out` with `hashNamesPhiValues.pl`, get only the geneID and phi mean columns that match those in `PopNames.out` and write them in csv format instead of tsv format.
    5. `phiValues009670-6.tsv`: Renamed from `009670-6.tsv` from Dr. Gilchrist's 2015 paper supplementary materials.
    6. `randGeneIDPhiMean.csv`: Use `rand500genes.pl` and `geneIDPhiMean.csv` to get 500 random genes that are in both Dr. Gilchrist's phi value file and Pop's RFP data file.
    7. `randNames.out`: Use `randGeneIDPhiMean.csv` as input into `printRandNames.pl` to print just the geneIDs of the random subset.
1. `codonTranslationRates.csv`: Raw codon translation rates from Pop's published work.
2. `orderedPopPADataRand500.csv`: `PopPADataRand500` ordered alphabetically to match easily with `orderedRandGeneIDPhiMean.csv`, produced by alphabetizeTwoFiles.pl`.
3. `orderedRandGeneIDPhiMean.csv`: `randGeneIDPhiMean.csv` ordered alphabetically to match easily with `orderedPopPADataRand500.csv`, produced by `alphabetizeTwoFiles.pl`.
4. `scripts`: Contains the scripts mentioned above to generate some preprocessedFiles.

