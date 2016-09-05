#!/usr/bin/env python2

import sys
import random
import glob
from rearrange_ms_dwgsim import get_refbases,get_altbases,get_strands
    
filelist = glob.glob(sys.argv[1]+'/*variants.txt')
filelist.append(sys.argv[1]+'/0.variants.txt')

#list of variant sites
with open(sys.argv[1]+'/allvars.txt') as f:
    #lines look like: chr22   39401   T       K
    varsites = [line.split()[1] for line in f]  #just pos

#invariant sites available to be hets       
availablesites = list(set(varsites) ^ set(range(1000000)))  #Return a new set with elements in either the set or other but not both.

#pick hets          
newvar = random.sample(availablesites, sys.argv[2])  #number of hets needed

#get ref/alt data for new hets
handle = open('ref/chr22_20-21M.fa', "rU")
ref = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()
refbases = get_refbases(ref['chr22'].seq,newvar)
refbases = [r.upper() for r in refbases]
altbases = get_altbases(refbases)  #returns list of hets for each variant position
strands = get_strands(sys.argv[2])

for f in filelist:
    #lines look like: chr22   39401   T       K       1
    f2 = open(f,'a')
    for i,b in enumerate(newvar):
        f2.write("chr22\t"+str(newvar)+"\t"+refbases[i]+"\t"+altbases[i]+"\t"+str(strands[i])+"\n")
    f2.close()
