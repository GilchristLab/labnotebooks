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
    print "Usage: python RandomGeneSelector.py csv_rfp_data csv_phi_values [NumGenes] [NumPositions]"
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
        numPositions = int(sys.argv[4])
    else:
        numPositions = -1
else:
    numGenes = -1
    numPositions = -1

#Total Codons and Genes
genes_phi = set()
genes_rfp = set()
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

data.close()
phis.close()

#Use all genes if no optional arguments
if numGenes == -1 or numGenes >= len(genes_rfp):
    numGenes = len(genes_rfp) - 1

#Randomly Sample Genes
genes = random.sample(genes_rfp, numGenes)

#New File names
if numPositions is -1: poses = "all"
else: poses = str(numPositions)
if numGenes is -1: gees = "all"
else: gees = str(numGenes)
rfpfile = "rfp_file_" + poses + "positions_" + gees + "genes.csv"
phifile = "phi_file_" + poses + "positions_" + gees + "genes.csv"

#Write RFP file
with open(rfpfile, "wb") as output:
    output.write(datatitle)
    rfpwriter = csv.writer(output,delimiter=",")
    for gene in genes:
        i = 0
        g = gs[gene]
        for pos in g.positions:
            if i == numPositions:
                break
            rfpwriter.writerow([gene, pos.pos, pos.codon, pos.rfpcount])
            i += 1

#Write Phi file
with open(phifile, "wb") as output:
    output.write(phistitle)
    phiwriter = csv.writer(output,delimiter=",")
    i = 1
    for gene in genes:
        g = gs[gene]
        phiwriter.writerow([gene, g.mean, g.median, g.ci_25, g.ci_975])

