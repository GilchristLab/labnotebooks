# input in format GeneID, Codon, Position, rfp_count
# Output should be in format GeneID, Position, Codon, RFPCount(s)

$datafile = "rfp.count.Etc.csv";

open(DATA, $datafile);

@genes = ();

while (chomp($line = <DATA>)){
    @data = split(/,/, $line);
# Trash header line
    if (@data[0] ne "GeneID"){
        $new_gene = {GeneID => $data[0], Codon => $data[1], 
                    Position => $data[2], rfpCount => $data[3]};
        push(@genes, $new_gene);
    }
}

print "GeneID,Position,Codon,rfpCount\n";

for $gene (@genes){
    print "$$gene{'GeneID'},$$gene{'Position'},$$gene{'Codon'},$$gene{'rfpCount'}\n";
}

close $datafile;
