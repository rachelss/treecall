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

N = len(vcffile.samples)
alignment = {s:[] for s in vcffile.samples}

for v in vcffile:
    if v.REF in bases and v.ALT[0] in bases:
        for s in vcffile.samples:
            geno = list(v.genotype(s).gt_bases)[::2] #genotype gets call for record
            alignment[s].append(''.join(geno))

output_handle = open(sys.argv[1]+'.alignment.phylip', "w")
alignment = [SeqRecord(Seq(''.join(alignment[s]), generic_dna), id=s) for s in sorted(alignment.keys())]
#print(alignment)
SeqIO.write(alignment, output_handle, "phylip")
output_handle.close()