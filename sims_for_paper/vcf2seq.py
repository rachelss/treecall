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
Usage: %s <input> [format]
<input> : a vcf or vcf.gz or '-' for stdin
[format]: fasta/phylip, default fasta
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
        out_fmt = sys.argv[2]
    else:
        out_fmt = 'fasta'

    vcffile.open()

    N = len(vcffile.samples)
    line = vcffile.fmt_line
    gt_codes = []
    while True:
        if line[0] == '#':
            pass
        else:
            v = Vcf(line)
            if v.REF == 'N' or v.ALT == '' or v.extract_info('INDEL') is not None:
                pass
            else:
                gt_codes.append(get_gt(v, vcffile))
        try:
            line = vcffile.next()
        except StopIteration:
            break
    vcffile.close()

    seqs = map(''.join, np.array(gt_codes).transpose().tolist())
    names = vcffile.samples
    seqiter = ('\t'.join([names[i],seqs[i]]) for i in xrange(N))

    from Bio import AlignIO
    alignments = AlignIO.parse(seqiter, 'tab')
    AlignIO.write(alignments, sys.stdout, out_fmt)
