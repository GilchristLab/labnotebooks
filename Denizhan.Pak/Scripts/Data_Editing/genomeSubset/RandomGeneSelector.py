import sys
import csv
import random

class Gene:
    def __init__(self, name):
        self.id = name
        self.codons = []

class Codon:
    def __init__(self, rfp_count, num_codons, codon):
        self.rfp_count = rfp_count
        self.num_codons = num_codons
        self.codon = codon

#Usage
if len(sys.argv) < 2:
    print ("Usage: python RandomGeneSelector.py csv_rfp_data [NumGenes]")
    exit(0)

#Open Arguments
data = open(sys.argv[1], 'r')

#Turn files in csvs
datatitle = data.readline()
dataReader = csv.reader(data, delimiter = ",")

#Check for optional arguments
if len(sys.argv) > 2:
    numGenes = int(sys.argv[2])
else:
    numGenes = -1

#Total Codons and Genes
genes_rfp = set()
genes_by_name = dict()

#Adds Genes from rfp file to rfp dataset
for row in dataReader:
    if len(row) != 4: continue
    name = row[0]
    if name not in genes_rfp:
        g = Gene(name)
        genes_rfp.add(name)
        genes_by_name[name] = g
    else:
        g = genes_by_name[name]
    rfp_count = row[1]
    codon_count = row[2]
    codon_name = row[3]
    g.codons.append(Codon(rfp_count, codon_count, codon_name))

data.close()

#Use all genes if no optional arguments
if numGenes == -1 or numGenes >= len(genes_rfp):
    numGenes = len(genes_rfp)

#Randomly Sample Genes
gene_names = random.sample(genes_rfp, numGenes)

#New File names
rfpfile = "simulated_rfp_file_" + str(numGenes) + "_genes.csv"

#Write RFP file
with open(rfpfile, "w") as output:
    output.write(datatitle)
    rfpwriter = csv.writer(output,delimiter=",")
    for gene_name in gene_names:
        gene = genes_by_name[gene_name]
        for codon in gene.codons:
            rfpwriter.writerow([gene_name, codon.rfp_count, codon.num_codons, codon.codon])
