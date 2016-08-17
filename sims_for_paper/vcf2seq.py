#!/usr/bin/env python2

import sys
import vcf
from Bio import AlignIO,SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

########################

bases = ['A','C','G','T']

#vcf has to have . not letters in info
vcffile = vcf.Reader(open(sys.argv[1], 'r'))

alignment = {s:[] for s in vcffile.samples}
alignment_vars = {s:[] for s in vcffile.samples}

for v in vcffile:
    if v.REF in bases and v.ALT[0] in bases:
        for s in vcffile.samples:
            geno = list(v.genotype(s).gt_bases)[::2] #genotype gets call for record; gets with a / in middle; use ::2 to remove
            alignment[s].extend(''.join(geno))
           
for i,j in enumerate(alignment[vcffile.samples[0]]):  #go through each site
    b = [alignment[s][i] for s in vcffile.samples]  #base for each sample at this site
    if len(set(b))>1: #check if variable
        for s in vcffile.samples:
            alignment_vars[s].append(alignment[s][i])  #put variable site into new seq for each sample
            
output_handle = open(sys.argv[1]+'.alignment.phylip', "w")
alignment_vars_list = [SeqRecord(Seq(''.join(alignment_vars[s]), generic_dna), id=s) for s in vcffile.samples]
SeqIO.write(alignment_vars_list, output_handle, "phylip")
output_handle.close()