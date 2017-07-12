
#input should be format GeneID, Position, Codon, rfp_count, from PopPAData.csv
#output will be to PopNames.out

$datafile = "../PopPAData.csv";
$outfile = "../PopNames.out";

open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

$lastGene = "";

@genes = ();

while (chomp($line = <DATA>)){
    @data = split(/,/, $line);
# Trash header line
    if (@data[0] ne "GeneID"){
# First line: Initialize all
        if ($lastGene eq ""){
        }

# New gene: Add name only
        elsif ($data[0] ne $lastGene){
            $new_gene = {ORF => $lastGene}; 
            push(@genes, $new_gene); 
        }

# Same gene: ignore
        $lastGene = $data[0];
    }
}

$new_gene = {ORF => $lastGene}; 
push(@genes, $new_gene);

print OUT "GeneID\n";

for $gene (@genes){
    print OUT "$$gene{'ORF'}\n";
}

close $datafile;
close $outfile;

