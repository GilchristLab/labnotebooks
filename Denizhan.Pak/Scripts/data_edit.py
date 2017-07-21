import sys
import csv

if sys.argv[1] == "-t":
    with open(sys.argv[2], "rb") as table:
        cr = csv.reader(table, dialect = csv.excel_tab)
	commaout = csv.writer(open('sample.csv', 'wb'), dialect=csv.excel)
	for line in cr:
	    commaout.writerow(line)
    exit()
else:
    genes = set()
    to_print = []
    fields = []

    with open(sys.argv[2], "rb") as phiSet:
        phiValues = csv.DictReader(phiSet)
        for row in phiValues:
            genes.add(row["ORF"])

    with open(sys.argv[1], "rb") as dataSet:
        data = csv.DictReader(dataSet)

        for row in data:
            if row["GeneID"] in genes:
                to_print.append(row)

    with open("EditedData.csv", "wb") as output:
        fieldnames = ["GeneID","Position","Codon","rfpCount"]
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        for line in to_print:
            writer.writerow(line)
