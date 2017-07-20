import sys
import csv
import random

class Gene:
    def __init__(self, name, mean, median, ci25, ci975):
        self.id = name
        self.codons = []
        self.positions = []
        self.mean = mean
        self.median = median
        self.ci_25 = ci25
        self.ci_975 = ci975

class Position:
    def __init__(self, numpos, cod, gen, rfp):
        self.pos = numpos
        self.codon = cod
        self.gene = gen
        self.rfpcount = rfp

#Usage
if len(sys.argv) < 3:
    print "Usage: python Genome_Data_Editor.py csv_rfp_dat csv_phi_values [NumGenes] [NumCodons]"
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
    if len(sys.argv) > 4:
        numCodons = int(sys.argv[4])
    else:
        numCodons = -1
else:
    numeGenes = -1
    numCodons = -1

#Total Codons and Genes
genes_phi = set()
genes_rfp = set()
codons = set()
gs = dict()

#Adds Genes with phi values to phi_set
for row in phisReader:
    if row[0] not in genes_phi:
        g = Gene(row[0], row[1], row[2], row[3], row[4])
        genes_phi.add(row[0])
        gs[row[0]] = g

#Add to rfp set if set phi set defined
for row in dataReader:
    if row[0] in genes_phi:
        g = gs[row[0]]
        genes_rfp.add(row[0])
        g.positions.append(Position(row[1], row[2], row[0],row[3]))

    #Add Codons
    if row[2] not in codons:
        codons.add(row[2])

data.close()
phis.close()

#Use all codons and genes if no optional arguments
if numCodons == -1 or numCodons >= len(codons):
    numCodons = len(codons) - 1
if numGenes == -1 or numGenes >= len(genes_rfp):
    numGenes = len(genes_rfp) - 1

#Randomly Sample Genes and codons
genes = random.sample(genes_rfp, numGenes)
codons = random.sample(codons, numCodons)

#New File names
rfpfile = "rfp_file_" + str(numCodons) + "codons_" + str(numGenes)+ "genes.csv"
phifile = "phi_file_" + str(numCodons) + "codons_" + str(numGenes)+ "genes.csv"

#Write RFP file
with open(rfpfile, "wb") as output:
    output.write(datatitle)
    rfpwriter = csv.writer(output,delimiter=",")
    for gene in genes:
        g = gs[gene]
        for pos in g.positions:
            if pos.codon in codons:
                rfpwriter.writerow([gene, pos.pos, pos.codon, pos.rfpcount])

#Write RFP file
with open(phifile, "wb") as output:
    output.write(phistitle)
    phiwriter = csv.writer(output,delimiter=",")
    i = 1
    for gene in genes:
        g = gs[gene]
        phiwriter.writerow([gene, g.mean, g.median, g.ci_25, g.ci_975])

