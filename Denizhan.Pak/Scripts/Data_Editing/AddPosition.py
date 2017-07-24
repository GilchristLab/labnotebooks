import sys
import csv
import random

class Gene:
    def __init__(self, name, mean, median, ci25, ci975):
        self.id = name
        self.codons = []
        self.positions = []

class Position:
    def __init__(self, numpos, cod, gen, rfp):
        self.pos = numpos
        self.codon = cod
        self.gene = gen
        self.rfpcount = rfp

with open(sys.argv[1]) as data:
    genes = set()
    sequence = []

