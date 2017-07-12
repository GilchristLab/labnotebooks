
# Script sorts the file below by Codon rather than its current random state (keeping gene the same)
$datafile = "SimulatedRFPData.csv";

open(DATA, $datafile);

@genes = ();

while (chomp($line = <DATA>)){
    @data = split(/,/, $line);
# Trash header line
    if (@data[0] ne "ORF"){
        $new_gene = {ORF => $data[0], rfpCount => $data[1], 
                    codonCount => $data[2], codon => $data[3]};
        push(@genes, $new_gene); 
    }
}

@sorted_genes = sort { 
    if ($$a{"ORF"} eq $$b{"ORF"}) {
        return $$a{"codon"} cmp $$b{"codon"};
    }            
} @genes;

for $gene (@sorted_genes){
    print "$$gene{'ORF'},$$gene{'rfpCount'},$$gene{'codonCount'},$$gene{'codon'}\n";
}

close $datafile;

