import sys

if len(sys.argv) != 3:
    print "Usage: CodonValCheck [name_of_data_file] [testing_output_file]"
    sys.exit()

codons = {}

with open(sys.argv[1], "r") as rfpData:
    rfpData.readline()
    for line in rfpData:
        codon = line[0:3]
        val = round(float(line[4:]), 4)
        codons[codon] = val

with open(sys.argv[2], "r") as outData:
    for line in outData:
        codon = line[0:3]
        val = round(float(line[4:]), 4)
        if codon in codons and codons[codon] != val:
            print codon + " is not Accurrate"
        if codon not in codons:
            print codon

print "All good"
