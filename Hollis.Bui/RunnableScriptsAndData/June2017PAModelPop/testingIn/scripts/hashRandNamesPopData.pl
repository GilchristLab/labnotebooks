
# inputs are 1) randNames.out and 2) PopPAData.csv
# Script will make a hash table comparing randNames.out with Pop's data's geneIDs
# saving the line from PopPAData.csv if the geneID matches.
# output should be to PopPADataRand500.csv, in format GeneID,Position,Codon,rfpCount

$namefile = "../preprocessedFiles/randNames.out";
$datafile = "../preprocessedFiles/PopPAData.csv";
$outfile = "../preprocessedFiles/PopPADataRand500.csv";

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

