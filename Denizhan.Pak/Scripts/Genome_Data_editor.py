import sys
import csv
import random

class Gene:
    def __init__(self, name):
        self.id = name
        self.codons = []
        self.positions = []
class Position:
    def __init__(self, numpos, cod, gen, rfp):
        self.pos = numpos
        self.codon = cod
        self.gene = gen
        self.rfpcount = rfp

if len(sys.argv) < 3:
    print "Usage: python Genome_Data_Editor.py csv_rfp_dat csv_phi_values [NumGenes] [NumCodons]"
    exit(0)

data = open(sys.argv[1], 'rb')
phis = open(sys.argv[2], 'rb')

if len(sys.argv) > 3:
    numGenes = int(sys.argv[2])
    if len(sys.argv) > 4:
        numCodons = int(sys.argv[3])
    else:
        numCodons = -1
else:
    numeGenes = -1
    numCodons = -1

dataReader = csv.reader(data, delimiter = ",")
phiReader = csv.reader(data, delimiter = ",")
genes = set()
codons = set()
l = []
gs = dict()
i = 0


for row in dataReader:
    if row[0] not in genes:
        g = Gene(row[0])
        genes.add(row[0])
        gs[row[0]] = g
        l.append(g)
    if row[2] not in codons:
        codons.add(row[2])
    g.positions.append(Position(row[1], row[2], row[0],row[3]))

if numCodons == -1 or numCodons >= len(codons):
    numCodons = len(codons) - 1
if numGenes == -1 or numGenes >= len(genes):
    numGenes = len(genes) - 1

genes = random.sample(gs.values(), numGenes)
codons = random.sample(codons, numCodons)

for gene in genes:
    i = 0
    for pos in gene.positions:
        i += 1
        if pos.codon in codons:
            print ", ".join([gene.id, pos.pos, pos.codon, pos.rfpcount])
