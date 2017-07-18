
# inputs are 1) PopPADataRand500.csv and 2) randGeneIDPhiMean.csv
# Script will read in both files then print them back out in alphabetical order.
# output should be to 1) orderedPopPADataRand500.csv and 2) orderedRandGeneIDPhiMean.csv

$PopFile = "../preprocessedFiles/PopPADataRand500.csv";
$PopOutFile = "../orderedPopPADataRand500.csv";
open(POP, $PopFile);
open(OUTPOP, ">", $PopOutFile); # Open for output

@genes = ();

while (chomp($line = <POP>)){
    push(@genes, $line);
}

@sorted_genes = sort(@genes);

for $gene (@sorted_genes){
    print OUTPOP "$gene\n";
}

close $PopFile;
close $PopOutFile;

##############################

$PhiMeanFile = "../preprocessedFiles/randGeneIDPhiMean.csv";
$PhiMeanOutFile = "../orderedRandGeneIDPhiMean.csv";
open(PHIMEAN, $PhiMeanFile);
open(OUTPHIMEAN, ">", $PhiMeanOutFile); # Open for output

@genes2 = ();

while (chomp($line = <PHIMEAN>)){
    push(@genes2, $line);
}

@sorted_genes = sort(@genes2);

for $gene (@sorted_genes){
    print OUTPHIMEAN "$gene\n";
}

close $PhiMeanFile;
close $PhiMeanOutFile;
