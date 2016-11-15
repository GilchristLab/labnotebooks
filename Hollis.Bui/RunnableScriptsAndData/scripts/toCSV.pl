# This script takes in the Lareau data and prints it back out in csv format.
# To modify the source file, change the name here.
$datafile = "../../../gilchrist-lfs/lareau/GSM1406453_untreated_1.percodon.txt";

# Output is to stdout, use redirection (>) to a file.

open(DATA, $datafile);

foreach $line (<DATA>){
    @data = split(/\s+/, $line);

    foreach $print (@data){
    # Print in CSV format
        if ($print eq $data[0]){
            print "$print";
        }
        else{
            print ",$print";
        }
    }

    print "\n";
}

close $datafile;
