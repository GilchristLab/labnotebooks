import sys
import csv
import random

class Gene:
    def __init__(self, name):
        self.id = name
        self.codons = dict()
        self.positions = []
        self.rfp_counts = []
        self.total = 0

class Codon:
    def __init__(self, cod, count, rfp):
        self.codon = cod
        self.count = int(count)
        self.rfpcount = int(rfp)
        self.positions = []

if len(sys.argv) != 2:
    print "usage python AddPosition.py rfp_file"
    exit()

with open(sys.argv[1]) as fdata:
    genes = dict()
    codons = set()
    fdata.readline()
    data = csv.reader(fdata, delimiter = ",")

    for line in data:
        name = line[0]
        if name not in genes:
            genes[name] = Gene(name)
        g = genes[name]
        g.codons[line[3]] = Codon(line[3], line[2], line[1])
        codons.add(line[3])
        g.total += int(line[2])

for g in genes:
    gene = genes[g]
    for i in range(gene.total):
        c = random.choice(list(codons))
        codon = gene.codons[c]
        while c.count < 1:
            c = random.choice(list(codons))
            codon = gene.codons[c]
        codon.count -= 1
        codon.positions.append(i)
        gene.positions.append(codon)
        gene.rfp_counts.append(0)

for g in genes:
    gene = genes[g]
    for name in codons:
        codon = gene.codons[name]
        if len(codon.positions) == 1:
            gene.rfp_counts[codon.positions[0]] = codon.rfpcount
        elif len(codon.positions) > 0:
            while codon.rfpcount > 0:
                index = codon.positions[random.randint(0, len(codon.positions) - 1)]
                rfp = random.randint(0, codon.rfpcount)
                gene.rfp_counts[index] += rfp
                codon.rfpcount -= rfp

if len(sys.argv) > 2:
    fname = sys.argv[2]
else:
    fname = "simRFPdata.csv"

with open(fname, "w") as output:
    output.write("GeneID,Position,Codon,rfpCount\n")
    for g in genes:
        i = 0
        gene = genes[g]
        for pos in gene.positions:
            output.write(gene.id + ", " + str(i) + ", " + str(pos.codon) + ", " + str(gene.rfp_counts[i]) + "\n")
            i += 1
