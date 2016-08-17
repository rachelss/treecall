#!/usr/bin/env python2

from __future__ import print_function
import sys
from os.path import basename as basename
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL) 
import vcf

def usage():
    print('''
Usage: %s <input> <outfile>
<input> : a vcf or vcf.gz
''' % basename(sys.argv[0]))
    sys.exit(0)

iupac_lookup = {
    'AA':'A', 'CC':'C','GG':'G','TT':'T',
    'AG':'R', 'CT':'Y', 'CG':'S',
    'AT':'W', 'GT':'K', 'AC':'M',
    'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V'
}

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
    
    outfile = open(sys.argv[2],'w')  
        
    bases = ['A','C','G','T']
    vcffile = vcf.Reader(open(sys.argv[1], 'r'))
    for v in vcffile:
        if v.REF in bases and v.ALT[0] in bases:  #check snp
            s = [str(b) for b in v.ALT if str(b) in bases] #filter X - ie list of alt bases
            s.insert(0,str(v.REF)) #put ref in front of base list
            
            #dict to convert gt numbers to bases
            if len(s) == 2:
                find_geno = {0:s[0]+s[0], 1:''.join(sorted(s[0]+s[1])), 2:s[1]+s[1]}
            elif len(s) == 3:
                find_geno = {0:s[0]+s[0], 1:''.join(sorted(s[0]+s[1])), 2:s[1]+s[1], 3:''.join(sorted(s[0]+s[2])), 4:''.join(sorted(s[1]+s[2])), 5:s[2]+s[2]}
            elif len(s) == 4:
                find_geno = {0:s[0]+s[0], 1:''.join(sorted(s[0]+s[1])), 2:s[1]+s[1], 3:''.join(sorted(s[0]+s[2])), 4:''.join(sorted(s[1]+s[2])), 5:s[2]+s[2], 6:''.join(sorted(s[0]+s[3])), 7:''.join(sorted(s[1]+s[3])), 8:''.join(sorted(s[2]+s[3])), 9:s[3]+s[3]}

            gt = [iupac_lookup[find_geno[v.genotype(sample).gt_type]] for sample in vcffile.samples]  #this variant; each sample; iupac of genotype as 0/1/2 converted to bases
            print(gt)    
              
            outfile.write('\t'.join((v.CHROM,str(v.POS),v.REF,'\t'.join(gt))))
    outfile.close()

    vcffile.close()
