
# inputs are 1) PopNames.out and 2) phiValues009670-6.tsv
# Script will make a hash table comparing PopNames.out with the phi values' names, 
# saving the phi mean if the name matches.
# output should be to geneIDPhiMean.csv, in format geneID,phiMean

$namefile = "../PopNames.out";
$datafile = "../phiValues009670-6.tsv";
$outfile = "../geneIDPhiMean.csv";

open(NAMES, $namefile);
open(DATA, $datafile);
open(OUT, ">", $outfile); # Open for output

%smallerPhiList = ();

while (chomp($line = <NAMES>)){
    $smallerPhiList{$line}{name} = $line;
}

print OUT "geneID,phiMean\n";

foreach $line (<DATA>){
    @data = split(/\t/, $line);
    $id = $data[0];

    if (exists($smallerPhiList{$id})){
        print OUT "$data[0],$data[1]\n";
    }
}

close $namefile;
close $datafile;
close $outfile;

