#!/usr/bin/env python2

import sys
import random
import glob
from Bio import AlignIO,SeqIO
from rearrange_ms_dwgsim import get_refbases,get_altbases,get_strands

#################################

def extract2(base):
    IUPAC = { 'M':('A','C'), 'R':('A','G'), 'W':('A','T'), 'S':('C','G'), 'Y':('C','T'), 'K':('G','T') }
    return(IUPAC[base])

#################################

probADO = float(sys.argv[2])
numhets = int(sys.argv[3])

#files for ado
inpath = sys.argv[1]+'HET/var'
filelist = glob.glob(inpath+'/*variants.txt')
filelist.remove(inpath+'/0.variants.txt')

for fname in filelist:
    f = open(fname,'r')
    f2 = open(fname.replace('HET','ADO'),'w')
    flines = f.readlines()
    #lines look like: chr22   39401   T       K       1
    #first lines aren't hets in root
    for line in flines[0:-numhets]:
        f2.write(line)
    #remaining lines are het in root - ADO allowed
    for line in flines[-numhets:]:
        #check if ADO
        p = random.random()
        print(p)
        if p<=probADO:
            splitline = line.split("\t")
            bps = extract2(splitline[3])
            splitline[3] = random.choice(bps) #pick one of two to keep (ie one to drop)
            f2.write("\t".join(splitline))
        else:
            f2.write(line)
    f2.close()
