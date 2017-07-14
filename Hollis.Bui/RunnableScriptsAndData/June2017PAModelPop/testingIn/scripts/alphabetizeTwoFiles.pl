
# inputs are 1) PopPADataRand500.csv and 2) randGeneIDPhiMean.csv
# Script will read in both files and then, using the same ordering as the first file,
# print out the matching phi mean lines in order to match up each gene with each other.
# output should be to orderedRandGeneIDPhiMean.csv, in format GeneID,phiMean

$PopFile = "../PopPADataRand500.csv";
$PhiMeanFile = "../randGeneIDPhiMean.csv";
$PopOutFile = "../orderedPopPADataRand500.csv";
$PhiMeanOutFile = "../orderedRandGeneIDPhiMean.csv";

open(POP, $PopFile);
open(PHIMEAN, $PhiMeanFile);
open(OUTPOP, ">", $PopOutFile); # Open for output
open(OUTPHIMEAN, ">", $PhiMeanOutFile); # Open for output

@genes = ();

while (chomp($line = <PHIMEAN>)){
    push(@genes, $line);
}

print OUTPHIMEAN sort(@genes);

@genes2 = ();

while (chomp($line = <POP>)){
    push(@genes2, $line);
}

print OUTPOP sort(@genes2);

close $PopFile;
close $PhiMeanFile;
close $PopOutFile;
close $PhiMeanOutFile;
