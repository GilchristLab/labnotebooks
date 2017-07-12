import sys
import csv

totalRFP = 0

class Gene:
    def __init__(self, name):
        self.id = name
        self.codons = []
        self.positions = []
        self.totRFP = int(0)

class Position:
    def __init__(self, numpos, cod, gen, rfp):
        global totalRFP
        self.pos = numpos
        self.codon = cod
        self.gene = gen
        self.rfpcount = rfp
        if rfp > 0:
            totalRFP += rfp
            gen.totRFP += rfp
        self.probRFP = 0.0


data = open(sys.argv[1], 'rb')

dataReader = csv.reader(data, delimiter = ",")
genes = set()
codons = set()

l = []
gs = dict()

for row in dataReader:
    if row[0] not in genes:
        g = Gene(row[0])
        genes.add(row[0])
        gs[row[0]] = g
        l.append(g)
    if row[2] not in codons:
        codons.add(row[2])
    g.positions.append(Position(row[1], row[2], gs[row[0]], int(row[3])))

print "Gene id, position, RFPprob"
for gene in gs.values():
    for pos in gene.positions:
        pos.probRFP = float(pos.rfpcount) / float(totalRFP)
        print gene.id + " " + str(pos.pos) + " " + str(pos.probRFP)
