import sys,Bio.Alphabet,Bio.Alphabet.IUPAC as IUPAC
from digest import *
from Bio.Seq import Seq

ins = open('proteinIn.txt').read().split()

peptide = Seq(ins[0],IUPAC.protein)
if not Bio.Alphabet._verify_alphabet(peptide):
    print 'Bad protein'
    exit(0)
    
enzyme = ins[6]
misses = ins[5]
minlen = ins[1]
maxlen = ins[2]
minweight = ins[3]
maxweight = ins[4]

if enzyme.upper() == 'CHYMOTRYPSIN':
    digester = Chymotrypsin(peptide,misses,minlen,maxlen,minweight,maxweight)
elif enzyme.upper() == 'CHYMOTRYPSIN_LOW':
    digester = Chymotrypsin_low(peptide,misses,minlen,maxlen,minweight,maxweight)
elif enzyme.upper() == 'TRYPSIN':
    digester = Trypsin(peptide,misses,minlen,maxlen,minweight,maxweight)

peptides = digester.digest()


####################### Output for CMD ##########################################
print 'Input seq:\t',digester.peptide
print 'Input Len:\t',len(peptide)
print 'Enzyme:\t\t',enzyme
print 'Max Missed:\t',misses
print 'Min length:\t',minlen
print 'Max length:\t',maxlen
print 'Min Weight:\t',minweight
print 'Max Weight:\t',maxweight
print

print 'Cleave Sites:\t', len(digester.sites)
print 'Positions:',
for s in digester.sites:
    print s+1,
print
print

for pep in sorted(peptides,key=lambda x: x[3]):
    print 'Peptide:\t\t', pep[0]
    print 'Length:\t\t\t', pep[1]
    print 'Molecular Weight:\t',pep[2]
    print 'Missed Cleaves:\t\t',pep[3]
    print 'Amino Left:\t\t',pep[4]
    print 'Amino Right:\t\t',pep[5]
    print









