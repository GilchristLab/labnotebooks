import sys
import csv
import random

class Gene:
    def __init__(self, name, phi):
        self.id = name
        self.codons = []
        self.positions = []
        self.phi = phi

class Position:
    def __init__(self, numpos, codon, geneID, rfp):
        self.pos = numpos
        self.codon = codon
        self.gene = geneID
        self.rfpcount = rfp

#Usage
if len(sys.argv) < 3:
    print ("Usage: python RandomGeneSelector.py csv_rfp_data csv_phi_values [NumGenes] [NumPositions]")
    exit(0)

#Open Arguments
data = open(sys.argv[1], 'r')
phis = open(sys.argv[2], 'r')

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
name_phi = dict()
genes = dict()

#Adds Genes with phi values to phi_set
for row in phisReader:
    geneID = row[0]
    phi = float(row[1])
    if geneID not in name_phi:
        name_phi[geneID] = phi

#Add to rfp set if set phi set defined
for row in dataReader:
    geneID = row[0]
    
    if geneID in genes:
        gene = genes[geneID]
        gene.positions.append(Position(row[1], row[2], row[0],row[3]))
    
    elif geneID in name_phi:
        gene = Gene(geneID, name_phi[geneID])
        gene.positions.append(Position(row[1], row[2], row[0],row[3]))
        genes[geneID] = gene

        
data.close()
phis.close()

#Use all genes if no optional arguments
if numGenes == -1 or numGenes >= len(name_phi):
    numGenes = len(genes) - 1

#Randomly Sample Genes
sampled_genes = random.sample(genes.keys(), numGenes)

#New File names
if numPositions is -1: poses = "all"
else: poses = str(numPositions)
if numGenes is -1: gees = "all"
else: gees = str(numGenes)
rfpfile = "rfp_file_" + poses + "positions_" + gees + "genes.csv"
phifile = "phi_file_" + poses + "positions_" + gees + "genes.csv"

#Write RFP file
with open(rfpfile, "w") as output:
    output.write(datatitle)
    rfpwriter = csv.writer(output,delimiter=",")
    for geneID in sampled_genes:
        i = 0
        gene = genes[geneID]
        for pos in gene.positions:
            if i == numPositions:
                break
            rfpwriter.writerow([geneID, pos.pos, pos.codon, pos.rfpcount])
            i += 1

#Write Phi file
with open(phifile, "w") as output:
    output.write(phistitle)
    phiwriter = csv.writer(output,delimiter=",")
    i = 1
    for geneID in sampled_genes:
        gene = genes[geneID]
        phiwriter.writerow([geneID, gene.phi])

