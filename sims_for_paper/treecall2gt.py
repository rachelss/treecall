#!/usr/bin/env python2
#treecall output looks like
#chr22  56101   G   1.88e-01    7.86e-01    GG          1.88e-01    GG          GT      6.28e-01    15      s3,s8,s2,s1
#chrom  pos     ref null_p      mut_p       null_base   null_base_p mut_base    mut_alt mut_conf_p  mut_loc mut_smpl

from __future__ import print_function
import sys

iupac_lookup = {
    'AA':'A', 'CC':'C', 'GG':'G', 'TT':'T',
    'AG':'R', 'CT':'Y', 'CG':'S',
    'AT':'W', 'GT':'K', 'AC':'M',
    'CGT':'B', 'AGT':'D', 'ACT':'H', 'ACG':'V'
}

def usage():
    print('''
usage: {} <input> <output> <n_sample>
<input>: treecall output, or '-' for stdin
'''.format(sys.argv[0]))
    return


if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()
        sys.exit(0)

    n_sample = int(sys.argv[3])

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1])
    
    outfile = open(sys.argv[2],'w')

    for line in f:
        chrom,pos,ref,null_p,mut_p,null_base,null_base_p,mut_base,mut_alt,mut_conf_p,mut_loc,mut_smpl = line.rstrip().split()
        base = iupac_lookup[mut_base]
        alt = iupac_lookup[mut_alt]
        mut_sidx = map(int, mut_smpl.replace('s','').split(','))
        gt = [base]*n_sample
        for i in mut_sidx:
            gt[i] = alt
            
        outfile.write('\t'.join((chrom,pos,base,'\t'.join(gt))))

    f.close()
    outfile.close() 
