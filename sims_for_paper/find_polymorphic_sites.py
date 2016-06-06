#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import vcf
import numpy as np
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def usage(mesg=None):
        if mesg is not None:
            print(mesg)
        print('''
Usage: %s <input> [filter]

<input> : a vcf or vcf.gz file or '-' for stdin
[filter]: AD:2;PL:80 # min_AD=2 && min_PL=80, default
          AD:5       # min_AD=5
          AD4:1       # min_AD4=1
          PL:90       # min_PL=90
''' % os.path.basename(sys.argv[0]))
        sys.exit(0)

def parse_filter_str(filter_str):
    filters = {}
    try:
        for f in filter_str.split(';'):
            k,v = f.split(':')
            filters[k] = int(v)
    except:
        usage('incorrect [filter] format')
    return filters

def pass_filter_by_AD(vcf_line, min_AD, samples):
    """
    FOR A SINGLE VARIANT
    Get AD (formerly DPR (="Number of high-quality bases observed for each allele"))
    Sum AD across samples for ref/alt(1)
    Check for het (ref and alt are present) for each sample - F=het; T=homo
    
    Return:
        Bool: T if both alleles have more than min reads AND both het and homo present
    """
    ad=[vcf_line.genotype(s).data.AD for s in samples]  #for each sample get allele depths
    ad=np.array(ad).astype(np.int)  #convert to np array of int
    
    s = np.sum(ad, axis=0)[0:2] # total dp for REF and top ALT across all samples
    het = np.logical_xor(ad[:,0], ad[:,1]) # presence of both alleles or not
    
    return np.min(s) >= min_AD and np.unique(het).size == 2

def pass_filter_by_AD4(vcf_line, min_AD4, samples):
    """
    AD4 = AD at ref F, alt F, ref R, alt R
    
    check at least F or R has min for both refs and alts
    
    """
    adf = [vcf_line.genotype(s).data.ADF for s in samples]  #for each sample get allele depth forward reads
    adr = [vcf_line.genotype(s).data.ADR for s in samples]  #for each sample get allele depth reverse reads
    adf = np.array(adf).astype(np.int)  #convert to np array of int
    adr = np.array(adr).astype(np.int)  #convert to np array of int

    goodf = adf >= min_AD4  #true/false array of >min [[T T] [T T]]
    goodr = adr >= min_AD4 
    ref_good = np.logical_and(goodf[:,0], goodr[:,0]) #check ref >min for all samples [T T]
    alt_good = np.logical_and(goodf[:,1], goodf[:,1]) #check alt>min for all samples [T T]
    
    return ref_good.sum()>0 and alt_good.sum()>0 #check at least F or R has min for both refs and alts

def pass_filter_by_PL(vcf_line, min_PL, samples):
    """
    Check if sum of PLs of most common genotype is > min required
    if not, most common geno is too common and not many hets
    if so, there's a reasonable prob that there's more than one geno across all samples
    
    Args:
        vcf_line: a variant
        min_PL: min phred scaled likelihood (summed across samples) required to pass filter
        samples: list of sample names
        
    Returns:
        bool: whether variant passes filter based on PL
    """
    pl=[vcf_line.genotype(s).data.PL for s in samples]  #for each sample get genotypes' pl
    pl=np.array(pl).astype(np.int)  #convert pl to np array of int
    
    most_common_geno_col = np.argmin(np.sum(pl>0, axis=0)) #get column of most common genotype by summing num times geno didn't occur
    penalty = np.sum(pl[:,most_common_geno_col]) #get sum of PLs of most_common_geno
    
    return penalty >= min_PL

def pass_filter(vcf_line, filters, samples):
    """
    Given requested functions and associated parameter value, check if vcf line passes
    """
    func = {'AD':pass_filter_by_AD, 'AD4':pass_filter_by_AD4, 'PL':pass_filter_by_PL}
    p = True
    for filt_name,filt_value in filters.iteritems():
        p = p and func[filt_name](vcf_line, filt_value, samples)
    return p


if __name__ == '__main__':

    if len(sys.argv) < 2:
        usage()
        
    if len(sys.argv) > 2:
        filters = parse_filter_str(sys.argv[2]) or {'AD':5,'PL':80}
    else:
        filters = {'AD':5,'PL':80}
        
    input_vcf = sys.argv[1]

    vcffile = vcf.Reader(open(input_vcf, 'r'))
    vcf_writer = vcf.Writer(open(input_vcf+'.vcf', 'w'), vcffile)

    samp_num = 0
    for v in vcffile:
        samp_num += 1
        if samp_num % 1000 == 0:
                print(str(samp_num) + ' ')
        if v.REF == 'N' or v.ALT == '' or v.is_indel:
            pass
        else:
            if pass_filter(v, filters,vcffile.samples):
                vcf_writer.write_record(v)