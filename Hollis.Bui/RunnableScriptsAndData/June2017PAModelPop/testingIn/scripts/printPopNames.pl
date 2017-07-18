
#input should be format GeneID,Position,Codon,rfpcount, from PopPAData.csv
#output will be to PopNames.out

$datafile = "../preprocessedFiles/PopPAData.csv";
$outfile = "../preprocessedFiles/PopNames.out";

open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

$lastGene = "";

@genes = ();

while (chomp($line = <DATA>)){
    @data = split(/,/, $line);
# Trash header line
    if (@data[0] ne "GeneID"){
        if ($lastGene eq ""){
# Skip to $lastGene = $data[0];
        }

# New gene: Add name only
        elsif ($data[0] ne $lastGene){
            $new_gene = {geneID => $lastGene}; 
            push(@genes, $new_gene); 
        }

# Same gene: ignore
        $lastGene = $data[0];
    }
}

$new_gene = {geneID => $lastGene}; 
push(@genes, $new_gene);

print OUT "GeneID\n";

for $gene (@genes){
    print OUT "$$gene{'geneID'}\n";
}

close $datafile;
close $outfile;

