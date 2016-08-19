#!/usr/bin/env python2

#output file is lines of
#chr22   140     G       G       R       G       G       G       G       G       G       G       G

from __future__ import print_function
import sys
import numpy as np

def usage():
    print('''
usage: {} <path> <n_sample> <n_snpsite>
'''.format(sys.argv[0]))

def init_matrix(dwgsim_input, n_sample, n_snpsite):
    """
    make matrix to track genotypes for samples
    
    args:
        file with lines of
            chr22   4901    T       K       
            chr     pos     ref     alt_het 
        number of samples (int)
        number of snps (int)
    
    return:
        dict: site/list of ref allele (len = n_sample+1 to include root sample)
    """
    true_gtype = {}
    with open(dwgsim_input) as f:
        for line in f:
            chrom,pos,ref,alt = line.rstrip().split()
            pos = int(pos)
            true_gtype[(chrom,pos)] = [ref]*(n_sample+1)
    return true_gtype

if __name__ == '__main__':
    if len(sys.argv) < 4:
        usage()
        sys.exit(0)

    path = sys.argv[1]
    n_sample = int(sys.argv[2])
    n_snpsite = int(sys.argv[3])

    allvars = path + '/var/allvars.txt'

    gt = init_matrix(allvars, n_sample, n_snpsite)
    for i in range(n_sample):
        sample_input = path+'/var/'+str(i+1)+'.variants.txt'
        with open(sample_input) as f:
            for line in f:
                chrom,pos,ref,alt,code = line.rstrip().split()
                pos = int(pos)
                assert (chrom,pos) in gt, 'site not found'
                gt[(chrom,pos)][i+1] = alt

    outfile = open(path + '/var/true.spgt.txt','w')
    for k in sorted(gt):
        linetowrite = [k[0], str(k[1])]
        linetowrite.extend(gt[k])
        outfile.write('\t'.join(linetowrite))
        outfile.write("\n")
    outfile.close()
