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

DELTA=0.0001  #move this so it's not global

def init_tree(tree):
    """
    node.sid = list of children

    """
    tree.leaf_order = map(int, tree.get_leaf_names())

    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            node.sid = [int(node.name)]
        else:
            node.name = ''
            node.sid = []
            for child in node.children:
                node.sid.extend(child.sid)

    m = len(tree)
    for i,node in zip(xrange(2*m-1), tree.traverse(strategy='postorder')):
        node.nid = i
        node.sid = sorted(node.sid)
        
    return tree

def p2phred(x):
    return -10.0*np.log10(x)

def phred2p(x):
    return 10.0**(-x/10.0)

def sum_PL(x, axis=None):
    return p2phred(phred2p(x).sum(axis=axis))

def normalize_PL(x):
    p = 10.0**(-x/10.0)
    return -10.0*np.log10(p/p.sum())

def normalize2d_PL(x):
    """
    Args:
        np.array (longdouble): PLs for a sample for all vars
        
    Returns:
         np.array (longdouble): PLs for a sample for all vars - rescaled slightly based on sum of all ll
    """
    p = 10.0**(-x/10.0)
    r = -10.0*np.log10(p/p.sum(axis=1)[:,None])
    return r

def gtype_distance(gt):
    """
    Args:
        gt(np.array (str)): genotypes as 1d array - usually either GTYPE3 (generic het/homos) or GTYPE10 (all possible gtypes)
    
    Return:
        np.array(int): Levenshtein (string) distance between pairs - eg AA-RR = 2
    """ 
    n = len(gt)
    gt_dist = np.zeros((n,n), dtype=int)
    for i,gi in enumerate(gt):
        for j,gj in enumerate(gt):
            gt_dist[i,j] = min(int(strdist(gi,gj)),int(strdist(gi,gj[::-1])))
            
    return gt_dist

def make_mut_matrix(mu, gtypes):
    """Makes a matrix for genotypes - only depends on mu
    
    Args:
        mu (int): mutation rate in Phred scale, default 80
        gtypes(np.array (str)): genotypes as 1d array - usually either GTYPE3 (generic het/homos) or GTYPE10 (all possible gtypes)
        
    Returns:
        np.array(float): substitution rate matrix
        np.array(float): substitution rate matrix with non-diagonal set to 0
        np.array(float): substitution rate matrix with diagonal set to 0
    """
    pmu = phred2p(mu)  #80 -> 10e-08
    gt_dist = gtype_distance(gtypes) #np.array: Levenshtein (string) distance between pairs - eg AA-RR = 2
    mm = pmu**gt_dist
    np.fill_diagonal(mm, 2.0-mm.sum(axis=0))
    mm0 = np.diagflat(mm.diagonal()) # substitution rate matrix with non-diagonal set to 0
    mm1 = mm - mm0 # substitution rate matrix with diagonal set to 0
    
    return mm,mm0,mm1

def make_mut_matrix_gtype3(mu):
    """same as above assuming gtype3 and w correct string distance for double mutation"""
    
    pmu = phred2p(mu)
    nmu = 1-pmu

    mm = np.array([[nmu**2, (2*pmu*nmu), (pmu*pmu)],
              [(nmu*pmu), (pmu**2)+(nmu**2), (nmu*pmu)], 
              [pmu*pmu, 2*pmu*nmu, (1-pmu)**2]])    
    
    mm0 = np.diagflat(mm.diagonal()) # substitution rate matrix with non-diagonal set to 0
    mm1 = mm - mm0 # substitution rate matrix with diagonal set to 0
    
    return mm,mm0,mm1

def make_mut_matrix_gtype10(mu):
    """same as above assuming gtype3 and w correct string distance for double mutation"""
    
    pmu = phred2p(mu)
    nmu = 1-pmu

    ##FIX ME
    mm = np.array([[nmu**2, (2*pmu*nmu), (pmu*pmu)],
              [(nmu*pmu), (pmu**2)+(nmu**2), (nmu*pmu)], 
              [pmu*pmu, 2*pmu*nmu, (1-pmu)**2]])    
    
    mm0 = np.diagflat(mm.diagonal()) # substitution rate matrix with non-diagonal set to 0
    mm1 = mm - mm0 # substitution rate matrix with diagonal set to 0
    
    return mm,mm0,mm1

def make_base_prior(het, gtypes):
    """Base prior probs
    for het=30, GTYPE3 = np.array(('RR','RA','AA'))
        [ 3.0124709,  33.012471,  3.0124709]

    for het=30, GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
        [ 6.0271094, 36.027109, 36.027109, 36.027109, 6.0271094, 36.027109, 36.027109, 6.0271094, 36.027109, 6.0271094]
    
    Args:
        het (int): heterozygous rate in Phred scale, default 30
        gtypes(np.array (str)): genotypes as 1d array
    
    Returns:
        np.array
    
    """
    return normalize_PL(np.array([g[0]!=g[1] for g in gtypes], dtype=np.longdouble)*het)

def calc_mut_likelihoods(tree, mm0, mm1):
    """
    go through tree from leaves to root - attach PLm to each node (not tips!)
    
    Args:
        tree (Tree)
        mm0: mutation matrix (np array of float) (non-diagonal set to 0)
        mm1: mutation matrix (np array of float) (diagonal set to 0)
        
    Returns:
        Tree (w annotated nodes)
    """
    n,g = tree.PL0.shape  #n = num var; g = num genos (eg 3)
    for node in tree.traverse(strategy='postorder'):
        if not node.is_leaf():
            node.PLm = np.zeros((2*len(node)-2,n,g), dtype=np.longdouble)  #len(node) = num tips associate

    for node in tree.traverse(strategy='postorder'):
        i = 0
        for child in node.children:
            sister = child.get_sisters()[0]
            if not child.is_leaf():
                l = child.PLm.shape[0]
                node.PLm[i:(i+l)] = p2phred(np.dot(phred2p(child.PLm), mm0)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
                i += l
            node.PLm[i] = p2phred(np.dot(phred2p(child.PL0), mm1)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
            i += 1

    return tree

def populate_tree_PL(tree, PLs, mm, attr): #e.g. populate_tree_PL(tree, PLs, mm0, 'PL0')
    """
    
    Args:
        tree (Tree)
        PLs (np.array): phred scaled likelihoods
        mm: mutation matrix (np array of float) (mm0 has non-diagonal set to 0; mm1 has diagonal set to 0)
        attr: attribute to be set e.g. PL0
    
    Returns:
        Tree: now has matrix attached to nodes
            PLs for all vars for this leaf or dot product of child's matrix and mutation matrix
    """
    n,m,g = PLs.shape # n sites, m samples, g gtypes
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            setattr(node, attr, PLs[:,node.sid[0],])  #sid is list of children's labels (numbers) - using 0 b/c only one label for leaf
        else:
            setattr(node, attr, np.zeros((n,g), dtype=np.longdouble))
            for child in node.children:
                setattr(node, attr, getattr(node, attr) + p2phred(np.dot(phred2p(getattr(child, attr)), mm))) #sum of phred of each child's likelihoods*mut matrix
                
    return tree