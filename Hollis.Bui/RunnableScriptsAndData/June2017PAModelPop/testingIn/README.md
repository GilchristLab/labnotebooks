# testing Input files

Contents:

0. `PopNames.out`: `PopPAData.csv` processed through `printPopNames.pl`: Just the geneIDs of the file.
1. `PopPAData.csv`: Processed version of the Pop data with RFPcounts (aka rfp.count... from Dr. Gilchrist's MATLAB code.
2. `codonTranslationRates.csv`: Raw codon translation rates from Pop's published work.
3. `geneIDPhiMean.csv`: Just the Phi mean values from `phiValues009670-6.tsv` that apply to `PopNames.out`, in .csv format.
4. `phiValues009670-6.tsv`: Originally `009670-6.tsv` from Dr. Gilchrist's 2015 paper supplementary materials.
5. `scripts`: Contains the scripts mentioned above to generate some of these files:
    0. `hashNamesPhiValues.pl`: Using the names from `PopNames.out` with the phi (mean) values from `phiValues009670-6.tsv1`, create `geneIDPhiMean.csv`.
    1. `printPopNames.pl`: Perl script to create `PopNames.out` from `PopPAData.csv`.

