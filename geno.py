#!/usr/bin/env python2

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>
# Author: Rachel Schwartz <Rachel dot Schwartz at asu dot edu>
# Author: Kael Dai <Kael dot Dai at asu dot edu>

from __future__ import print_function
import warnings
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import sys
import itertools
import numpy as np
from scipy.stats import sem
from editdistance import eval as strdist
import vcf

with warnings.catch_warnings(ImportWarning):
    from ete2 import Tree

warnings.filterwarnings('error')

from tree_est import *

from utils import *

def read_vcf_records(vcffile, fmt, maxn=1000):
    """Read vcf file - get info about variants - need to clarify how this is different from read_vcf

    Args:
        vcffile (VcfFile)
        fmt(dict): item:pos_in_list from the 9th column of a vcf line
        maxn (int): number of lines / sites in file to process

    Returns:
        np.array (tuple): variant info (chrom, pos, ref)
        np.array (int): Number of high-quality bases observed for each of the 2 most common alleles
        np.array (double): List of Phred-scaled genotype likelihoods for each of the 2 most common alleles

    """    
    print('read next %d sites' % maxn, end=' ', file=sys.stderr)
    variants,DPRs,PLs = [],[],[]
    i = 0
    for v in vcffile:
        i += 1
        try:
            variants.append((v.CHROM,v.POS,v.REF))
            dpr = np.array(v.extract_gtype('DPR', fmt, v.get_DPR4), dtype=np.uint16)
            DPRs.append(dpr)
            pl = np.array(v.extract_gtype('PL', fmt, v.get_PL10), dtype=np.longdouble)
            PLs.append(pl)
        except Exception as e:
            print(e, file=sys.stderr)
            print(v, file=sys.stderr)
        if i == maxn:
            print('... %s:%s ...' % (v.CHROM, v.POS), end=' ', file=sys.stderr)
            break

    variants = np.array(variants)
    DPRs = np.array(DPRs, dtype=np.uint16)
    PLs = np.array(PLs, dtype=np.longdouble)

    print(' done', file=sys.stderr)
    return variants, DPRs, PLs

def genotype_main(args):
    """
    uses init_tree, make_base_prior, make_mut_matrix, read_vcf_records, genotype
    """
    
    GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
    print(args, file=sys.stderr)

    tree = Tree(args.tree)
    tree = init_tree(tree)

    base_prior = make_base_prior(args.het, GTYPE10) # base genotype prior
    mm,mm0,mm1 = make_mut_matrix(args.mu, GTYPE10) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    vcffile = VcfFile(args.vcf)
    fmt = vcffile.fmt
    fout = open(args.output, 'w')
    fout.close()
    fout = open(args.output, 'a')
    score = 0
    while True:
        variants, DPRs, PLs = read_vcf_records(vcffile, vcffile.fmt, args.nsite)
        records,s = genotype(PLs, tree, variants, mm, mm0, mm1, base_prior)
        np.savetxt(fout, records, fmt=['%s','%d','%s','%.2e','%.2e','%s','%.2e','%s','%s','%.2e','%d','%s'], delimiter='\t')
        score += s
        if len(PLs) < args.nsite:
            break
    print('sum(PL) = %.2f' % score)
    fout.close()

def genotype(PLs, tree, variants, mm, mm0, mm1, base_prior):
    """
    uses populate_tree_PL, calc_mut_likelihoods, phred2p
    """
    GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
    # calculate total likelihoods for each genotypes
    populate_tree_PL(tree, PLs, mm, 'PL') # dim(tree.PL) = site x gtype
    tree_PL = tree.PL + base_prior
    # calculate no-mutation likelihoods for each genotypes
    #try:
    populate_tree_PL(tree, PLs, mm0, 'PL0') # dim(tree.PL0) = site x gtype
    #except Exception as e:
    #    print('populate_tree_PL():', e, file=sys.stderr)
    #    sys.exit(1)
    tree_PL0 = tree.PL0 + base_prior
    # calculate mutation likelihoods for each genotypes and mutation locations
    calc_mut_likelihoods(tree, mm0, mm1)
    mut_PLs = np.swapaxes(tree.PLm,0,1) # site x location x gtype
    mut_PLs += base_prior
    n,l,g = mut_PLs.shape # n sites, l locations, g gtypes
    nn = np.arange(n)

    k = tree_PL.argmin(axis=1) # most likely base genotype for each site
    tree_P_per_site = phred2p(tree_PL).sum(axis=1) # total tree likelihood

    k0 = tree_PL0.argmin(axis=1) # most likely non-mutation base genotype for each site
    null_PL = tree_PL0[nn,k0] # best non-mutation likelihood (across genotypes) for each site
    null_P_per_site = phred2p(tree_PL0).sum(axis=1) # total non-mutation likelihood

    k1 = np.array([np.unravel_index(s.argmin(), (l,g)) for s in mut_PLs]) # site x 2, most likely mutation event for each site
    k1l = k1[:,0] # most likely location
    k1g = k1[:,1] # most likely base genotype
    mut_PL = mut_PLs[nn,k1l,k1g] # best mutation likelihood (across location and genotypes) for each site
    mut_P_per_site = phred2p(mut_PLs).sum(axis=(1,2)) # total mutation likelihood

    null_PLs = np.array([node.PL0 for node in tree.iter_descendants(strategy='postorder')])
    k2 = null_PLs[k1l,nn,].argmin(axis=-1) # get most likely mutation mutant genotype

    node_sids = np.array([','.join(map(str,node.sid)) for node in tree.iter_descendants(strategy='postorder')])
    records = np.array(zip(
            variants[nn,0],                         # chrom
            variants[nn,1],                         # pos
            variants[nn,2],                         # ref
            null_P_per_site/tree_P_per_site,        # null_P
            mut_P_per_site/tree_P_per_site,         # mut_P
           #GTYPE10[k],                             # MLE_base_gtype
           #phred2p(tree_PL[nn,k])/tree_P_per_site, # MLE_base_gtype_P
            GTYPE10[k0],                            # MLE_null_base_gtype
            phred2p(null_PL)/tree_P_per_site,       # MLE_null_base_gtype_P
            GTYPE10[k1g],                           # MLE_mut_base_gtype
            GTYPE10[k2],                            # MLE_mut_alt_gtype
            phred2p(mut_PL)/tree_P_per_site,        # MLE_mut_base_gtype_P
            k1l,                                    # MLE_mut_location
            node_sids[k1l]),                        # MLE_mut_samples
        dtype=[
            ('chrom','a10'),('pos','i4'),('ref','a1'),
            ('null_p','f8'),('mut_p','f8'),
            ('null_base','a2'),('null_base_p','f8'),
            ('mut_base','a2'),('mut_alt','a2'),('mut_conf_p','f8'),
            ('mut_loc','i4'),('mut_smpl','a128')])
    score = p2phred(records['mut_p']+records['null_p']).sum()
    return records,score