#!/usr/bin/env python

from __future__ import print_function

import sys
from os.path import basename as basename
from string import maketrans
import numpy as np
from memoize import memoized
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL) 

from pyvcf import Vcf, VcfFile

def usage(mesg=None):
    if mesg is not None:
        print(mesg)
    print('''
Usage: %s <input> [sidx_range] [format]
<input> : a vcf or vcf.gz or '-' for stdin
[sidx_range]: sample index range, 0-based in pythonic style, e.g "0:4" means the first 3 samples
[format]: matrix/long, default matrix
''' % basename(sys.argv[0]))
    sys.exit(0)

fmt_choices = {'fasta', 'phylip'}
iupac_lookup = {
    'A,G':'R', 'C,T':'Y', 'C,G':'S',
    'A,T':'W', 'G,T':'K', 'A,C':'M',
    'C,G,T':'B', 'A,G,T':'D', 'A,C,T':'H', 'A,C,G':'V'
}
gt_delim = maketrans('|', '/')

@memoized
def make_gcode(ref, alt):
    global iupac_lookup
    allele = [ref] + alt.split(',')
    n = len(allele)
    gcode = []
    for i in xrange(n):
        for j in xrange(i+1):
            if i==j:
                gcode.append(allele[i])
            else:
                het = ','.join(sorted([allele[j],allele[i]]))
                gcode.append(iupac_lookup[het])
    return gcode


def get_gt(v, vf):
    global gt_delim
    if v.ALT == '':
        return [v.REF]*len(vf.samples)
    else:
        gcode = make_gcode(v.REF, v.ALT)
        if 'GT' in vf.fmt:
            gt = v.extract_gtype('GT', vf.fmt)
            idx = np.sum(np.array([g.translate(gt_delim).split('/') for g in gt]).astype(np.int), axis=1)
        elif 'PL' in vf.fmt:
            idx = np.argmin(np.array(v.extract_gtype('PL', vf.fmt, str.split, ',')).astype(np.int), axis=1)
        return [gcode[i] for i in idx]



if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage()

    vcffile = VcfFile(sys.argv[1])

    if len(sys.argv) > 2:
        sidx = sys.argv[2]
    else:
        sidx = None

    if len(sys.argv) > 3:
        out_fmt = sys.argv[3]
    else:
        out_fmt = 'matrix'

    vcffile.open()

    samples = vcffile.samples
    if sidx is not None:
        samples = eval('samples['+sidx+']',{'samples':samples,'os':None}) # precaution that 'os.system("rm -rf /")' won't happen
    #N = len(samples)
    sites = []
    gt_codes = []
    for v in vcffile:
        if v.REF == 'N' or v.ALT == '' or v.extract_info('INDEL') is not None:
            pass
        else:
            gt = get_gt(v, vcffile)
            if sidx is not None:
                gt = eval('gt['+sidx+']',{'gt':gt,'os':None}) # precaution that 'os.system("rm -rf /")' won't happen
            site = (v.CHROM,v.POS)
            if out_fmt == 'matrix':
                print('\t'.join((v.CHROM,v.POS,v.REF,'\t'.join(gt))))
            else:
                for s,g in zip(samples,gt):
                    print('\t'.join((v.CHROM,v.POS,v.REF,s,g)))
    vcffile.close()
