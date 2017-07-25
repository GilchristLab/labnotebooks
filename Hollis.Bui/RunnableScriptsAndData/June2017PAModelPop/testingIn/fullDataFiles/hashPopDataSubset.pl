
# inputs are 1) PopSubset.out and 2) PopPAData.csv
# Script will make a hash table comparing PopSubset.out with Pop's data's geneIDs
# saving the line from PopPAData.csv if the geneID matches.
# output should be to PopPADataSubset.csv, in format GeneID,Position,Codon,rfpCount

$namefile = "PopSubset.out";
$datafile = "../preprocessedFiles/PopPAData.csv";
$outfile = "PopPADataSubset.csv";

open(NAMES, $namefile);
open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

%smallerPopDataList = ();

while (chomp($line = <NAMES>)){
    $smallerPopDataList{$line}{name} = $line;
}

foreach $line (<DATA>){
    @data = split(/,/, $line);
    $id = $data[0];

    if (exists($smallerPopDataList{$id})){
        print OUT "$line";
    }
}

close $namefile;
close $datafile;
close $outfile;

