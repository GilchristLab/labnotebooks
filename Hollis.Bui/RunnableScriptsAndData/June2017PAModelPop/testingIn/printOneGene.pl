
#input should be from PopPADataSubset.csv
#output will be to PopPADataOneGene.csv 

$datafile = "fullDataFiles/PopPADataSubset.csv";
$outfile = "PopPADataOneGene.csv";

open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

$lastGene = "";

while (chomp($line = <DATA>)){
    @data = split(/,/, $line);
# Skip header line
    if (@data[0] ne "GeneID"){
# Print only the data that matches the first GeneID found
        if ($lastGene eq ""){
            $lastGene = $data[0];
        }

        if ($data[0] eq $lastGene){
            print OUT "$line\n";
        }
    }
# Print header line back out
    else{
        print OUT "$line\n";
    }
}

close $datafile;
close $outfile;

