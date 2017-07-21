# input should be geneIDPhiMean.csv, in format GeneID,phiMean.
# output should be to randGeneIDPhiMean.csv, in the same format but only as a subset.

$datafile = "../preprocessedFiles/geneIDPhiMean.csv";
$outfile = "../preprocessedFiles/randGeneIDPhiMean.csv";

open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

@genes = ();
while (chomp($line = <DATA>)){
    push(@genes, $line); 
}

@out = ();
for (1 .. 500){
    push(@out, splice @genes, rand @genes, 1);
}

print OUT "geneID,phiMean\n";

for $gene (@out){
    print OUT "$gene\n";
}

close $datafile;
close $outfile;

