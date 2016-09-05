#!/usr/bin/env python2

import sys
import random
import glob
from Bio import AlignIO,SeqIO
from rearrange_ms_dwgsim import get_refbases,get_altbases,get_strands
    
#files that new hets will be added to
path = sys.argv[1]
filelist = glob.glob(path+'/*variants.txt')
filelist.append(path+'/0.variants.txt')

#list of variant sites
with open(path+'/allvars.txt') as f:
    #lines look like: chr22   39401   T       K
    varsites = [line.split()[1] for line in f]  #just pos

#invariant sites available to be hets       
availablesites = list(set(varsites) ^ set(range(1000000)))  #Return a new set with elements in either the set or other but not both.

#pick hets
num_new_hets = int(sys.argv[2])
newvar = random.sample(availablesites, num_new_hets)  #number of hets needed

#get ref/alt data for new hets
handle = open('ref/chr22_20-21M.fa', "rU")
ref = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()
refbases = get_refbases(ref['chr22'].seq,newvar)
refbases = [r.upper() for r in refbases]
altbases = get_altbases(refbases)  #returns list of hets for each variant position
strands = get_strands(num_new_hets)

for f in filelist:
    #lines look like: chr22   39401   T       K       1
    f2 = open(f,'a')
    for i,b in enumerate(newvar):
        f2.write("chr22\t"+str(b+1)+"\t"+refbases[i]+"\t"+altbases[i]+"\t"+str(strands[i])+"\n")
    f2.close()
