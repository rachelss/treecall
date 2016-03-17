#!/usr/bin/env python

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>

from __future__ import print_function
import sys
from os.path import splitext
import gzip
from memoize import memoized

def warning(*obj):
    print(*obj, file=sys.stderr)


@memoized
def PL10_order(ref, alt):
    """Get the position of the genotype likelihood in the PL field

    There are 10 possible genotypes. The PL field in a vcf file has a normalized Phred-scaled likelihood for
    each of the possible genotypes. This provides the expected order of genotypes in that field.

    Args:
        ref (str): Reference base
        alt (str): Alternate base(s) (separated by ,)

    Returns:
        list: For each allele (all 10) its genotype number e.g. 0 for homozygous ref, 1 for het, 2 for homo alt

    """
    nt = ('A','C','G','T')
    assert ref in nt, "Reference allele not a nucleotide"
    
    if alt == '':
        a = [ref]
    else:
        alts = alt.split(',')
        for al in alts:
            assert al in nt, "Alternate allele not a nucleotide"
        a = [ref] + alt.split(',')  #a is a list of alleles
        
    aidx = {b:i for i,b in enumerate(a)}  #aidx is a dict of allele:pos in list

    for x in nt:
        if x not in aidx:
            aidx[x] = len(a)
    
    return [aidx[b1]+aidx[b2]*(aidx[b2]+1)/2 for i,b1 in enumerate(nt) for b2 in nt[i:]]


@memoized
def DPR4_order(ref, alt):
    """Get the position of the base in order of ref,alt,not

    Args:
        ref (str): Reference base
        alt (str): Alternate base(s) (separated by ,)

    Returns:
        list: 

    """
    nt = ('A','C','G','T')
    assert ref in nt, "Reference allele not a nucleotide"
    
    if alt == '':
        a = [ref]
    else:
        alts = alt.split(',')
        for al in alts:
            assert al in nt, "Alternate allele not a nucleotide"
        a = [ref] + alt.split(',')  #a is a list of alleles
            
    #[ unicode(x.strip()) if x is not None else '' for x in row ]
    return [a.index(x) if x in a else len(a) for x in nt]


class Vcf(object):
    #example vcf line:
    #1 898921 . C G,<X> 0 . DP=211;I16=112,51,1,0,7057,346989,16,256,8061,399183,50,2500,3220,73296,25,625;QS=0.997431,0.00256946,0;SGB=-0.379885;RPB=1;MQB=1;MQSB=0.792466;BQB=1;MQ0F=0    PL:DP:DV:DPR    0,255,255,255,255,255:164:1:163,1,0
    
    def __init__(self, line=None, fixed_only=False):
        if line is not None:
            self.line = line.rstrip('\n')
            self.lazy_parse(fixed_only)

    def lazy_parse(self, fixed_only=True):
        fields = self.line.split('\t')
        self.CHROM = fields[0]
        self.POS = fields[1]
        self.ID = fields[2]
        self.REF = fields[3]
        self.ALT = fields[4].rstrip('<X>').rstrip(',')
        self.QUAL = fields[5]
        self.FILTER = fields[6]
        self.INFO = fields[7]
        self.info = None
        if not fixed_only:
            self.FMT = fields[8]
            self.extra = fields[9:]
            self.gtypes = {}

    def extract_gtype(self, tag, fmt, func=None, *args):
        """Get info from col 10 of vcf line for particular criterion eg genotype likelihood or num bases
            add info for that criterion on self.gtypes
    
        Args:
            self (Vcf): a line of a vcf file with all the info about the variant
            tag (str): either 'DPR' or 'PL' (Number of high-quality bases observed for each allele or List of Phred-scaled genotype likelihoods)
            fmt (dict): item:pos_in_list from the 9th column of a vcf line e.g. PL:1 DP:2 DV:3 DPR:4
            func (func): apply this to the particular info from col 10 eg splits PL string by commas if func=str.split and args=','
            args (): see func
    
        Returns:
            list: values or adjusted values (eg genotype likelihoods) for particular tag from vcf line
    
        """
        if tag in self.gtypes:
            return self.gtypes[tag]
        gtype = []
        for s in self.extra:
            content = (s.split(':'))[fmt[tag]]  #get 10th col of vcf as list - access the values for item of interest e.g. PL or DPR
            if func is None:
                gtype.append(content)   #can use values
            else:
                gtype.append(func(content, *args))  #can adjust PL or DPR values - eg splits PL string by commas if func=str.split and args=','
        self.gtypes[tag] = gtype
        return gtype

    def extract_info(self, tag=None, func=None, *args):
        if self.info is None:
            info = {}
            for s in self.INFO.split(';'):
                k,eq,v = s.partition('=')
                if eq == '':
                    v = True
                else:
                    v = v.split(',')
                    if len(v) == 1:
                        v = v[0]
                info[k] = v
            self.info = info
        if tag is None:
            return self.info
        else:
            value = self.info.get(tag)
            if func is not None:
                return func(value, *args)
            else:
                return value

    def get_PL10(self, x):
        order = PL10_order(self.REF, self.ALT)
        pl = x.split(',')
        n = len(pl)
        return [pl[i] if i < n else '255' for i in order]

    def get_DPR4(self, x):
        order = DPR4_order(self.REF, self.ALT)
        dpr = x.split(',')
        if max(order) >= len(dpr):
            warning(self)
        return [dpr[i] for i in order]

    def __str__(self):
        """Return the vcf line"""
        return self.line


class VcfFile(object):
    def __init__(self, filename=None, ifmt=None):
        if filename is not None:
            self.filename = filename
            self.ifmt = ifmt or self._get_format(filename)      #either z for .gz or v (for vcf?)
            self.opened = False
            self.vcf_hdr = []
            self.samples = []
            self.fmt = {}   #a dict of item:pos_in_list from the 9th column of a vcf line
            self.read_vcf_header()
            if len(self.samples) > 0:
                self.read_FMT()     #sets fmt to a dict of item:pos_in_list

    def __iter__(self):
        self.open()
        if self.seekable:
            line = '#'
        else:
            line = self.fmt_line
        while True:
            if line[0] == '#':
                pass
            else:
                yield Vcf(line)
            try:
                line = self.next()
            except:
                self.close()
                break

    @staticmethod
    def parse_FMT(fmt):
        """Splits a string by : and converts it to a dict of item:pos in list  (applied to the 9th col of a vcf)"""
        tags = fmt.split(':')
        return {v:i for i,v in enumerate(tags)}

    def _get_format(self, filename=None):
        """Get format of vcf file - either z or v for .gz or not"""
        filename = filename or self.filename
        base, ext = splitext(filename)
        return 'z' if ext == '.gz' else 'v'

    def open(self, force=False):
        if self.opened:
            if not force:
                return self.f
            else:
                warning('re-opening %s' % self.filename)
                self.close()
        self.seekable = True
        if self.ifmt == 'z':
            f = gzip.open(self.filename, 'rb')
        else:
            if self.filename == '-':
                f = sys.stdin
                self.seekable = False
            else:
                f = open(self.filename, 'r')
        self.opened = True
        self.f = f

    def close(self):
        if self.opened:
            self.f.close()
            self.opened = False

    def next(self):
        return self.f.next()

    def read_vcf_header(self):
        self.open(force=True)
        for line in self.f:
            if line[0] == '#':
                self.vcf_hdr.append(line.rstrip('\n'))
                if line[1] != '#':
                    fields = line.rstrip('\n').split('\t')
                    self.samples = fields[9:]
                    if self.seekable:
                        self.data_pos = self.f.tell()
                    break
            else:
                raise ValueError('No header present.')

    def skip_header(self):
        if not self.opened:
            raise ValueError('I/O operation on closed file.')
        if self.seekable:
            self.f.seek(self.data_pos)
        else:
            pass

    def read_FMT(self):
        if not self.opened:
            self.open(force=True)
        while True:
            line = self.next()
            if line[0] != '#':
                break
        self.fmt_line = line.rstrip('\n')       #get first non-info line of vcf
        fields = self.fmt_line.split('\t')
        self.fmt = VcfFile.parse_FMT(fields[8]) #fields[8] is a bunch of info such as PL:DP:DV:DPR - parsing makes a dict of item:pos_in_list
        if self.seekable:
            self.f.seek(self.data_pos)
