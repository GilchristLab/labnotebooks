import sys
import csv
import random

class Gene:
    def __init__(self, name, phi):
        self.id = name
        self.codons = []
        self.phi = phi

class Codon:
    def __init__(self, rfp_count, num_codons, codon):
        self.rfp_count = rfp_count
        self.num_codons = num_codons
        self.codon = codon

#Usage
if len(sys.argv) < 3:
    print "Usage: python RandomGeneSelector.py csv_rfp_data csv_phi_values [NumGenes]"
    exit(0)

#Open Arguments
data = open(sys.argv[1], 'rb')
phis = open(sys.argv[2], 'rb')

#Turn files in csvs
datatitle = data.readline()
phistitle = phis.readline()
dataReader = csv.reader(data, delimiter = ",")
phisReader = csv.reader(phis, delimiter = ",")

#Check for optional arguments
if len(sys.argv) > 3:
    numGenes = int(sys.argv[3])
else:
    numGenes = -1

#Total Codons and Genes
genes_phi = set()
genes_rfp = set()
genes_by_name = dict()

#Adds Genes from phi file to phi_set
for row in phisReader:
    name = row[0]
    phi = row[1]
    if name not in genes_phi:
        g = Gene(name, phi)
        genes_phi.add(name)
        genes_by_name[name] = g

#Add genes to rfp file set if corresponding phi exists
for row in dataReader:
    name = row[0]
    if name in genes_phi:
        g = genes_by_name[name]

        rfp_count = row[1]
        codon_count = row[2]
        codon_name = row[3]

        genes_rfp.add(name)
        g.codons.append(Codon(rfp_count, codon_count, codon_name))

data.close()
phis.close()

#Use all genes if no optional arguments
if numGenes == -1 or numGenes >= len(genes_rfp):
    numGenes = len(genes_rfp)

#Randomly Sample Genes
gene_names = random.sample(genes_rfp, numGenes)

#New File names
rfpfile = "simulated_rfp_file_" + str(numGenes) + "_genes.csv"
phifile = "simulated_phi_file_" + str(numGenes) + "_genes.csv"

#Write RFP file
with open(rfpfile, "wb") as output:
    output.write(datatitle)
    rfpwriter = csv.writer(output,delimiter=",")
    for gene_name in gene_names:
        gene = genes_by_name[gene_name]
        for codon in gene.codons:
            rfpwriter.writerow([gene_name, codon.rfp_count, codon.num_codons, codon.codon])

#Write Phi file
with open(phifile, "wb") as output:
    output.write(phistitle)
    phiwriter = csv.writer(output,delimiter=",")
    for gene_name in gene_names:
        gene = genes_by_name[gene_name]
        phiwriter.writerow([gene_name, gene.phi])

