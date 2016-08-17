#!/usr/bin/env python

from __future__ import print_function
import sys
import itertools

def usage():
    print('''
usage: {} <matrix0> <matrix1> <n_snpsite>
<matrix0>: true gtypes
<n_snpsite>: number of true snp sites
format of <matrix>: chrom,pos,ref,gt1,gt2,...
'''.format(sys.argv[0]))
    return

def read_gt_matrix(filename, n_snp=None):
    """returns dict of site / list of gts """
    gtypes = dict();
    snps = dict();
    with open(filename) as f:
        for i,line in enumerate(f):
            fields = line.rstrip().split()
            chrom = fields.pop(0)
            pos = fields.pop(0)
            site = chrom+':'+pos
            gtypes[site] = fields 
            if n_snp is not None and i<n_snp:
                snps[site] = True
    return gtypes,snps


if __name__ == '__main__':
    if len(sys.argv) < 4:
        usage()
        sys.exit(0)

    n_snp = int(sys.argv[3])
    gtypes0,snps = read_gt_matrix(sys.argv[1], n_snp)
    gtypes1,empty = read_gt_matrix(sys.argv[2])
    n_sample = len(gtypes0.itervalues().next())-1

    site0 = set(gtypes0.keys())
    site1 = set(gtypes1.keys())
    n_site = len(site0)

    fp = 0  #false pos (mutation incorrectly estimated)
    fn = 0  #false neg (mutation not est)
    tp = 0  #true pos (mutations estimated)
    tn = 0  #true neg (no mutation, none found)
    for site in site0 - site1:  #true snps not est - some gt are correct (tn), some are not (fn)
        gt0 = gtypes0[site]
        ref = gt0.pop(0)
        if site in snps:
            x = sum([gt!=ref for gt in gt0])  #number of true mutations
            fn += x  #real mut but not found
            tn += n_sample - x  #not mut and none est
        else:
            tn += n_sample

    for site in site1 - site0:  #false pos est snps
        gt1 = gtypes1[site]
        ref = gt1.pop(0)
        x = sum([gt!=ref for gt in gt1])  #num that are dif from ref
        fp += x
        tn = n_sample - x

    for site in site0 & site1:  #true snps that were est
        gt0 = gtypes0[site]
        gt1 = gtypes1[site]
        ref0 = gt0.pop(0)
        ref1 = gt1.pop(0)
        if site in snps:
            for g0,g1 in zip(gt0,gt1):
                if g0==ref0:
                    if g1==g0:
                        tn += 1
                    else:
                        fp += 1
                else:
                    if g1==g0:
                        tp += 1
                    else:
                        fn += 1
        else:
            for g0,g1 in zip(gt0,gt1):
                if g1==g0:
                    tn += 1
                else:
                    fp += 1

    N = n_sample*n_site
    precision = float(tp)/(tp+fp)  # = positive predictive value
    recall = float(tp)/(tp+fn)  #recall = Sensitivity, true positive rate, or probability of detection
    print('{:.4f}\t{:.4f}'.format(precision,recall))  
