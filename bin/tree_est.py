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

def read_vcf(filename, evidence=60):
    """Read vcf file - get info about variants

    Args:
        filename (str): vcf filename
        evidence (int): minimum evidence in Phred scale
            for a site to be considered, default 60

    Returns:
        vcf: a vcffile
        np.array (tuple): variant info (chrom, pos, ref)
            for each variant
        np.array (int): allele depth for each of the 2 most common alleles
            for each variant overall; array for each sample for each variant
            in order of freq
        np.array (int): List of Phred-scaled genotype likelihoods
            for each of the genotypes from the 2 most common alleles for each variant
            by common/common, common/less_common, less/less

    """
    vcffile = vcf.Reader(open(filename, 'r'))
    bases = ['A','C','G','T']
    variants,ADs,PLs = [],[],[]
    
    #mapping allele to genotype eg row 0, col 2 is ref/alt2 = most common alleles
    #gives PLs 0,3,5
    #triallelic the PL pattern is RR,RA1,A1A1,RA2,A1A2,A2A2 
    #so correctly gets PLs for RR, RA2, A2A2
    a2g = np.array((
            ((0,0,0), (0,1,2), (0,3,5), (0,6,9)),
            ((2,1,0), (2,2,2), (2,4,5), (2,7,9)),
            ((5,3,0), (5,4,2), (5,5,5), (5,8,9)),
            ((9,6,0), (9,7,2), (9,8,5), (9,9,9))
        ))
    
    for v in vcffile:
        if v.REF in bases and v.ALT[0] in bases:
            variants.append((v.CHROM,v.POS,v.REF))
            
            #ad for each sample for each allele
            ad = np.array([v.genotype(s).data.AD for s in vcffile.samples], dtype=np.uint16) 
            #triallelic the PL pattern is RR,RA1,A1A1,RA2,A1A2,A2A2 
            pl = np.array([v.genotype(s).data.PL for s in vcffile.samples], dtype=np.uint16) #list of PL-as-int in array
            
            #ak = columns of two most common alleles ordered by freq
            #sum ad across samples; order by decreasing depth, take only the two most common alleles
            ak = ad.sum(axis=0).argsort(kind='mergesort')[-2:][::-1]   
            #append ad ordered by freq NOT ref,alt
            ADs.append(ad[...,ak])  
            
            #get genotypes' PLs in order of allele freq across samples
            gk = a2g[ak[0],ak[1]]
            PLs.append(pl[...,gk])
    
    variants = np.array(variants)
    ADs = np.array(ADs)
    PLs = np.array(PLs)
    
    #for each variant, sum PL for each genotype across samples
    #genotypes are ordered from most to least likely NOT ref/ref, ref/alt, alt/alt
    #check if PL sums are all greater than evidence
    #this removes sites where the joint genotyping likelihood across samples
    #for second most likely genotype is < 10^-6
    #i.e. most likely genotype for each sample has strong support
    #k_ev = (np.sort(PLs).sum(axis=1)>=evidence).sum(axis=1)==2  #this used to be ==3 but that doesn't seem right - it should be checked
    #variants,ADs,PLs = variants[k_ev],ADs[k_ev],PLs[k_ev]
    #commented about above bc it's probably filtering uncertain variants but that's the whole point of treecall, and we're not sure what it's doing anyway

    print(' done', file=sys.stderr)
    return vcffile, variants, ADs, PLs


def neighbor_main(args):
    """generate neighbor-joining tree then do recursive NNI and recursive reroot

    Args:
        vcf(str): input vcf/vcf.gz file, "-" for stdin
        output(str): output basename
        mu (int): mutation rate in Phred scale, default 80
            WHY IS THE DEFAULT 80????? IE 10^-8!!!!!!!!!!!!!!!!
        het (int): heterozygous rate in Phred scale, default 30
        min_ev(int): minimum evidence in Phred scale for a site to be considered
            default 60
    
    Output:
        newick trees
    
    """
   
    print(args, file=sys.stderr)
    vcffile, variants, DPRs, PLs = read_vcf(args.vcf, args.min_ev)
    #variants =  np.array (tuple): variant info (chrom, pos, ref)  for each variant
    #DPRs = np.array (int): Number of high-quality bases observed for each of the 2 most common alleles for each variant
    #PLs = np.array (int): List of Phred-scaled genotype likelihoods for each of the 2 most common alleles (3 genotypes) for each variant
    
    GTYPE3 = np.array(('RR','RA','AA'))
    base_prior = make_base_prior(args.het, GTYPE3) # base genotype prior; heterozygous rate in Phred scale, default 30; e.g. for het=30 [ 3.0124709,  33.012471,  3.0124709]
    mm,mm0,mm1 = make_mut_matrix_gtype3(args.mu) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)  # pairwise differences between samples based only on PLs (should include mutation, but also shouldn't matter)
    allscores = []
    for i in range(10):  #10 different starting trees
        tree = init_star_tree(n_smpl)
        internals = np.arange(n_smpl)
        
        #1st tree is nj tree (tho with raw scores not adjusted for saturation)
        if i==0:
            D,tree = neighbor_joining(D.copy(), tree.copy(), internals) #haven't checked this; make nj tree and update D given internal nodes; pass copy
        #all other trees are semi-random
        else:
            tree.set_outgroup(str(i))
            tree.resolve_polytomy()
    
        tree = init_tree(tree.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
        tree = populate_tree_PL(tree.copy(), PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
        tree = calc_mut_likelihoods(tree.copy(), mm0, mm1)  #add PLs w mutation
        
        tree.write(outfile=args.output+'.nj0.tre', format=5)
        
        rerooted = 1
        while rerooted > 0:
            best_tree,best_PL = recursive_NNI(tree.copy(), PLs, mm0, mm1, base_prior,DELTA)
            #print(best_tree)
            best_tree,best_PL,rerooted = recursive_reroot(best_tree.copy(), PLs,mm0, mm1, base_prior,DELTA)  #why are brlens negative?
            #print(best_tree)
            print('PL_per_site = %.4f' % (best_PL/n_site))
            best_tree.write(outfile=args.output+'.nj.'+str(i)+'.tre', format=5)
            allscores.append(best_PL)
        i+=1
    
    print(allscores)
    
def init_star_tree(n):
    """Creates a tree, adds n children in star with numbers as names

    Args:
        n (int): Number of children in tree

    Returns:
        Tree: 
    """
    
    tree = Tree()
    for i in xrange(n):
        tree.add_child(name=str(i))
    return tree


def pairwise_diff(PLs, i, j):
    """
    Calculates difference between pairs of samples - sum across variants - where dif for any given var depends on PL for each possible genotype
    For each var: Dif geno will be close to 1; Same geno will be close to 0
    But depends slightly on PL
    
    Returns:
        numpy.float128 (one value): sum of diffs across vars
    """
    
    pli = normalize2d_PL(PLs[:,i])  #adjust PLs for sample i slightly
    plj = normalize2d_PL(PLs[:,j])  #adjust PLs for sample j slightly
    p = phred2p(pli+plj) # n x g
    return (1-p.sum(axis=1)).sum()  


def make_D(PLs):
    """
    Get pairwise differences between samples based on PLs (e.g. for generating nj tree)
    
    Args:
        PLs (np.array (longdouble)): List of Phred-scaled genotype likelihoods
            for each of the 2 most common alleles for each variant
        
    Returns:
        np.array (longdouble): matrix of n x n where n = number nodes in bifurcating tree
            includes diff btwn tips but not internal nodes - computed during tree building
    """
    n,m,g = PLs.shape   #n_site,n_smpl,n_gtype
    D = np.zeros(shape=(2*m-2,2*m-2), dtype=np.longdouble) #2*m-2 = the total number of nodes of unrooted tree (incl leaves and internal)
    for i,j in itertools.combinations(xrange(m),2):
        D[i,j] = pairwise_diff(PLs, i, j)
        D[j,i] = D[i,j]
    return D


def neighbor_joining(D, tree, internals):
    #fsum will have better precision when adding distances across sites
    #based on PLs not mutation
    """
    
    Args:
        D (np.array): pairwise differences between samples based on PLs (passing copy)
        tree (Tree): tree of class Tree with num tips = num samples
        internals (np.array): array of sample numbers
        
    Returns:
        Tree
        D (np.array): update pairwise differences now there are internal nodes to compare
    
    """
    print('neighbor_joining() begin', end=' ', file=sys.stderr)
    m = len(internals)
    while m > 2:  #if m is 2 then only two connected to root
        d = D[internals[:,None],internals]  #initially D matrix w/o 0 distance btwn internal nodes; then add in nodes as they have distances
        u = d.sum(axis=1)/(m-2)

        Q = np.zeros(shape=(m,m), dtype=np.longdouble)
        for i,j in itertools.combinations(xrange(m),2):  #std Q matrix calc
            Q[i,j] = d[i,j]-u[i]-u[j]
            Q[j,i] = Q[i,j]
        #print(Q.astype(int))
        np.fill_diagonal(Q, np.inf)
        #print(np.unique(Q, return_counts=True))
        i,j = np.unravel_index(Q.argmin(), (m,m))  #location in matrix of smallest Q value (ie closest nodes/tips)
        l = len(D)+2-m

        for k in xrange(m):
            D[l,internals[k]] = D[internals[k],l] = d[i,k]+d[j,k]-d[i,j]
        D[l,internals[i]] = D[internals[i],l] = vi = (d[i,j]+u[i]-u[j])/2
        D[l,internals[j]] = D[internals[j],l] = vj = (d[i,j]+u[j]-u[i])/2

        ci = tree&str(internals[i])
        cj = tree&str(internals[j])
        ci.detach()
        cj.detach()
        node = Tree(name=str(l))
        node.add_child(ci,dist=int(vi))
        node.add_child(cj,dist=int(vj))
        tree.add_child(node)
        #print(tree)

        internals = np.delete(internals, [i,j])
        internals = np.append(internals, l)
        m = len(internals)
        print('.', end='', file=sys.stderr)

    print(' done', file=sys.stderr)
    return D,tree

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

def update_PL(node, mm0, mm1):
    """
    PL for nodes depend on children so must be updated if node children change due to nni/reroot
    
    node has
        PL0: PLs given no mutation
        PLm: w mutation
        children
        
    changes PLm and PL0 for node and all its children (recursive); also sid for children
    
    returns:
        Tree: could be subtree; includes all recursively updated PL0 and PLm and node labels (sid)
    """

    n,g = node.PL0.shape
    l = 2*len(node)-2
    #node.PL0 = np.zeros((n,g), dtype=np.longdouble)
    node.PL0.fill(0.0)
    node.PLm = np.zeros((l,n,g), dtype=np.longdouble)
    for child in node.children:
        sid = sorted(map(int,child.get_leaf_names()))
        if child.sid != sid:  #sid is supposed to be names of leaves - could be dif due to swapping in nni
            newchild = update_PL(child, mm0, mm1)
            newchild.sid = sid
            child.detach()
            node.add_child(newchild)
        node.PL0 += p2phred(np.dot(phred2p(child.PL0), mm0)) 
    i = 0
    for child in node.children:
        sister = child.get_sisters()[0]
        if not child.is_leaf():
            l = child.PLm.shape[0]
            node.PLm[i:(i+l)] = p2phred(np.dot(phred2p(child.PLm), mm0)) + p2phred(np.dot(phred2p(sister.PL0), mm0))
            i += l
        node.PLm[i] = p2phred(np.dot(phred2p(child.PL0), mm1)) + p2phred(np.dot(phred2p(sister.PL0), mm0)) 
        i += 1

    return node

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

def score(tree, base_prior):
    """
    used to compare trees
    
    Args:
        Tree
            PL0: PLs for all vars for leaf or dot product of child's matrix and mutation matrix
                where mm0: substitution rate matrix with non-diagonal set to 0
            PLm:
                    mm1: substitution rate matrix with diagonal set to 0
    
    """
    Pm = phred2p(tree.PLm+base_prior).sum(axis=(0,2))       #why add baseprior
    P0 = phred2p(tree.PL0+base_prior).sum(axis=1)
    return p2phred(Pm+P0).sum()


def annotate_nodes(tree, attr, values):

    for node in tree.iter_descendants('postorder'):
        setattr(node, attr, values[node.nid])

    return tree

def read_label(filename):
    """from tab delim file: dict of
        key: first col in file / index
        value: second col in file (or first col if only one)
    """
    label = {}
    with open(filename) as f:
        i = 0
        for line in f:
            c = line.rstrip().split('\t')
            if len(c) > 1:
                label[c[0]] = c[1]
            else:
                label[str(i)] = c[0]
            i += 1
    return label

def partition(PLs, tree, sidx, min_ev):
    if tree.is_root():
        print('partition() begin', end=' ', file=sys.stderr)
    m = len(sidx) # number of samples under current node
    print(m, end='.', file=sys.stderr)
    if m == 2:
        child1 = tree.add_child(name=str(sidx[0]))
        child1.add_features(samples=np.atleast_1d(sidx[0]))
        child2 = tree.add_child(name=str(sidx[1]))
        child2.add_features(samples=np.atleast_1d(sidx[1]))
    elif m > 2:
        smat = make_selection_matrix2(m)
        pt, cost = calc_minimum_pt_cost(PLs, smat, min_ev)
        k0 = pt==0
        sidx0 = np.atleast_1d(sidx[k0])
        child = tree.add_child(name=','.join(sidx0.astype(str)))
        child.add_features(samples=sidx0)
        if len(sidx0) > 1:
            partition(PLs[:,k0,], child, sidx0, min_ev)
        k1 = pt==1
        sidx1 = np.atleast_1d(sidx[k1])
        child = tree.add_child(name=','.join(sidx1.astype(str)))
        child.add_features(samples=sidx1)
        if len(sidx1) > 1:
            partition(PLs[:,k1,], child, sidx1, min_ev)
    else:
        print('m<=1: shouldn\'t reach here', file=sys.stderr)
        sys.exit(1)
    if tree.is_root():
        print(' done', file=sys.stderr)


def calc_minimum_pt_cost(PLs, smat, min_ev):
    n,m,g = PLs.shape
    pt_cost = np.inf
    for k in smat:
        x0 = PLs[:,k==0,].sum(axis=1) # dim = n_site x 2
        x0min = x0.min(axis=1) # dim = n_site x 1
        x0max = x0.max(axis=1) # dim = n_site x 1
        x1 = PLs[:,k==1,].sum(axis=1) # dim = n_site x 2
        x1min = x1.min(axis=1) # dim = n_site x 1
        x1max = x1.max(axis=1) # dim = n_site x 1
        # take everything
        #c = (x0 + x1).sum()
        # cap the penalty by mu
        #c = (x0>mu).sum()*mu + x0[x0<=mu].sum() + (x1>mu).sum()*mu + x1[x1<=mu].sum()
        # ignore sites where signal from either partitions is weak
        #c = (x0min+x1min)[(x0max>min_ev) & (x1max>min_ev)].sum()
        # ignore sites where signals from both partitions are weak
        c = (x0min+x1min)[(x0max>min_ev) | (x1max>min_ev)].sum()
        # some weird cost function that broadly penalize partition of similar samples
        #k0 = x0.argmin(axis=1)
        #k1 = x1.argmin(axis=1)
        #c = np.minimum(x0[k0],x1[k1]).sum() + (k0==k1).sum()*mu
        if c < pt_cost:
            pt_cost = c
            pt = k
    return pt, pt_cost


def make_selection_matrix(m, t=20):
    n = 2**(m-1)
    if m>3 and m<=t: # special treatment for intermediate size
        l = (m,)*n
        x = np.array(map(tuple, map(str.zfill, [b[2:] for b in map(bin, xrange(4))], (3,)*4)), dtype=np.byte)
        y = np.zeros((n,m),dtype=np.byte)
        for i in xrange(m-3):
            a,b = x.shape
            y[0:a,-b:] = x
            y[a:(2*a),-b:] = x
            y[a:(2*a),-b] = 1
            x = y[0:(2*a),-(b+1):]
        for s in y:
            yield s
    else:
        for i in xrange(n):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)


def make_selection_matrix2(m, t=20):
    n = 2**(m-1)
    if m>3 and m<=t: # special treatment for intermediate size
        l = (m,)*n
        x = np.array(map(tuple, map(str.zfill, [b[2:] for b in map(bin, xrange(4))], (3,)*4)), dtype=np.byte)
        y = np.zeros((n,m),dtype=np.byte)
        for i in xrange(m-3):
            a,b = x.shape
            y[0:a,-b:] = x
            y[a:(2*a),-b:] = x
            y[a:(2*a),-b] = 1
            x = y[0:(2*a),-(b+1):]
        for s in y:
            yield s
    elif m<=3:
        for i in xrange(n):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)
    else:
        r1 = np.random.randint(1,m-1,2**t)
        r2 = np.random.rand(2**t)
        x = ((1+r2)*2**r1).astype(int)
        for i in iter(x):
            yield np.array(tuple(bin(i)[2:].zfill(m)), dtype=np.byte)


def reroot(tree, PLs, mm0, mm1, base_prior,DELTA):
    """
    
    return:
        Tree
        np.array (PLs)
        int: flag if rerooted (1) or not (0)
    """
    '''
              /-A              /-A              /-B
           /-|              /-|              /-|
    -root-|   \-B => -root-|   \-C => -root-|   \-C
          |                |                |
           \-C              \-B              \-A
    '''

    best_tree = tree
    best_PL = score(tree, base_prior)
    flag = 0

    for node in tree.iter_descendants('postorder'):  #go through all nodes including tips but not root
        tree_reroot = tree.copy()
        new_root = tree_reroot.search_nodes(sid=node.sid)[0]  #gets node of interest
        tree_reroot.set_outgroup(new_root)  #sets node of interest to outgroup
        
        tree_reroot = init_tree(tree_reroot)  #tree has nid's (node id) and sid's (list of tip names - sorted)
        tree_reroot = populate_tree_PL(tree_reroot, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
        tree_reroot = calc_mut_likelihoods(tree_reroot, mm0, mm1)  #add PLs w mutation
        
        #tree_reroot = update_PL(tree_reroot, mm0, mm1)  #new PL given decendants
        PL_reroot = score(tree_reroot, base_prior) 
        #print(tree_reroot)
        #print(PL_reroot)
        if PL_reroot < best_PL * (1-DELTA): #new best tree only if significantly better ie trees could be similar but status quo wins
            best_tree = tree_reroot
            best_PL = PL_reroot
            flag = 1
            
    return best_tree,best_PL,flag


def recursive_reroot(tree, PLs,mm0, mm1, base_prior,DELTA):
    """
    starting at tips, work up tree, get best way of rooting subtree 
    """

    print('recursive_reroot() begin', file=sys.stderr)
    PL = score(tree, base_prior)
    for node in tree.iter_descendants('postorder'):  #go through all nodes including tips but not root
        rerooted = 0
        new_tree = tree.copy()
        node_leaves = node.get_leaf_names()
        new_node = new_tree.get_common_ancestor(node_leaves)  #get corresponding node in new tree
        try:
            print('Node:')
            print(node)
            new_tree.set_outgroup(new_node) #reroot
    
            #recalculate
            new_tree = init_tree(new_tree.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
            new_tree = populate_tree_PL(new_tree.copy(), PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
            new_tree = calc_mut_likelihoods(new_tree.copy(), mm0, mm1)  #add PLs w mutation
                   
            PL_new = score(new_tree, base_prior)
            
            print('This tree:')
            print(new_tree)
            print('Has this score:')
            print(PL_new)
    
            if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                best_tree = new_tree.copy()
                PL = PL_new
                rerooted = 1
        except:
            print('Can\'t set node as outgroup:')
            print(node)
            
    if rerooted == 1:  #there was a better tree
        tree = best_tree.copy()
        print('Best tree:')
        print(tree)
    else:                        
        print('No change to tree:')
        print(tree)
        print(PL)
        
    print(' done', end='', file=sys.stderr)
    print(tree)
    print(PL)
    return tree,PL,rerooted


def nearest_neighbor_interchange(node, PLs,mm0, mm1, base_prior,DELTA):
    '''
    Args:
    
    Return:
        node: Tree (could be a subtree of the full tree)
        np.array(PL)
        int: flag to indicate nni happened
        
    Process:
    
              /-A              /-A              /-A
           /-|              /-|              /-|
          |   \-B          |   \-C          |   \-D
    -node-|       => -node-|       => -node-|
          |   /-C          |   /-B          |   /-B
           \-|              \-|              \-|
              \-D              \-D              \-C
           ||               ||               ||
           \/               \/               \/
        reroot()         reroot()         reroot()
    '''
    
    c1,c2 = node.children  #children of root node
    possible_rearrangements = []
    
    #children are leaves - don't need to swap anything
    if c1.is_leaf() and c2.is_leaf():
        return [None]
    
    #one child is a leaf - rerooting will provide all possible combinations - flagged if rerooted
    elif c1.is_leaf():
        c21,c22 = c2.children
        node.set_outgroup(c22)
        possible_rearrangements.append(node.copy())
        node.set_outgroup(c21)
        possible_rearrangements.append(node.copy())
        return possible_rearrangements
        
    elif c2.is_leaf():
        c12,c11 = c1.children
        node.set_outgroup(c12)
        possible_rearrangements.append(node.copy())
        node.set_outgroup(c11)
        possible_rearrangements.append(node.copy())
        return possible_rearrangements

    else:
        #rerootings of original tree
        node_copy1 = node.copy()
        c1,c2 = node_copy1.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        for n in [c11,c12,c21,c22]:
            node_copy1.set_outgroup(n)
            possible_rearrangements.append(node_copy1.copy())

        #2nd tree - swap relationships and reroot
        node_copy2 = node.copy()
        c1,c2 = node_copy2.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        c12 = c12.detach()
        c22 = c22.detach()
        c1.add_child(c22)
        c2.add_child(c12)
        
        possible_rearrangements.append(node_copy2.copy())
        for n in [c11,c12,c21,c22]:
            node_copy2.set_outgroup(n)
            possible_rearrangements.append(node_copy2.copy())
            
        #3rd tree - swap relationships and reroot
        node_copy3 = node.copy()
        c1,c2 = node_copy3.children
        c11,c12 = c1.children
        c21,c22 = c2.children
        c12 = c12.detach()
        c21 = c21.detach()
        c1.add_child(c21)
        c2.add_child(c12)
        
        possible_rearrangements.append(node_copy3.copy())
        for n in [c11,c12,c21,c22]:
            node_copy3.set_outgroup(n)
            possible_rearrangements.append(node_copy3.copy())

        return possible_rearrangements


def recursive_NNI(tree, PLs, mm0, mm1, base_prior,DELTA):
    #recursive just means traverse the tree 
    """
    
    Args:
        tree(Tree)
        mm0: mutation matrix (np array of float) (non-diagonal set to 0)
        mm1: mutation matrix (np array of float) (diagonal set to 0)
        base_prior (np.array): Base prior probs depending on het pl

    Returns:
        Tree (tree)
        np.array (PL): phred-scaled likelihooods
        
    tree resulting from each round of nni (working from tips to root) printed to trees_tried.txt
    
    """
    print('recursive_NNI() begin', end=' ', file=sys.stderr)
    PL = score(tree, base_prior)
    #goes until can get through tree w/o nni at any node
    #a la phylip
    num_nnis=1
    print(tree)
    print(PL)
    while(num_nnis>0):
        num_nnis=0
        print('Start nni round')
#        for node in tree.traverse('postorder'):
#            print(node)
        for node in tree.traverse('postorder'):
            #goes through each node, does nni if better
            if node.is_leaf():
                continue
            print('.', end='', file=sys.stderr)
            print(node)
            possible_rearrangements = nearest_neighbor_interchange(node.copy(),PLs, mm0, mm1, base_prior,DELTA)
#            print('Original original tree:')
#            print(tree)
            if possible_rearrangements[0] is not None:
#                for r in possible_rearrangements:
#                    print(r)
                if node.is_root():  #because can't get parent as for below
                    for r in possible_rearrangements:
                        print('Rearranged node:')
                        print(r)
                        
                        new_tree = init_tree(r.copy())  #tree has nid's (node id) and sid's (list of tip names - sorted)
                        new_tree = populate_tree_PL(new_tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
                        new_tree = calc_mut_likelihoods(new_tree, mm0, mm1)  #add PLs w mutation

                        PL_new = score(new_tree, base_prior)
                        
                        print('This tree:')
                        print(new_tree)
                        print('Has this score:')
                        print(PL_new)

                        if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                            best_tree = new_tree.copy()
                            print(best_tree)
                            PL = PL_new
                            num_nnis = 1
                    if num_nnis == 1:  #there was a better tree
                        tree = best_tree.copy()
                        print('Best tree so far:')
                        print(tree)
                    else:                        
                        print('No change to tree:')
                        print(tree)
                        print(PL)

                else:
                    for r in possible_rearrangements:
                        print('Rearranged node:')
                        print(r)
                        new_tree = tree.copy()
                        node_leaves = node.get_leaf_names()
                        new_node = new_tree.get_common_ancestor(node_leaves)  #get corresponding node in new tree
                        
                        parent = new_node.up
                        new_node.detach()
                        parent.add_child(r)
                        
                        #modify copy of tree
                        new_tree = init_tree(new_tree)  #tree has nid's (node id) and sid's (list of tip names - sorted)
                        new_tree = populate_tree_PL(new_tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
                        new_tree = calc_mut_likelihoods(new_tree, mm0, mm1)  #add PLs w mutation
                        
                        PL_new = score(new_tree, base_prior)
                                                
                        print('This tree:')
                        print(new_tree)
                        print('Has this score:')
                        print(PL_new)
                        
                        if PL_new < (PL-DELTA): #should this be multiplied or subtracted?
                            PL = PL_new
                            num_nnis = 1
                            best_tree = new_tree.copy()
#                        print('Original tree:')
#                        print(tree)

                    if num_nnis == 1:  #there was a better tree
                        tree = best_tree.copy()
                        print('Best tree so far:')
                        print(tree)
                        break  #take best tree and start over because now nni's will be all different
                    else:                        
                        print('No change to tree:')
                        print(tree)
                        print(PL)
                    
        #print(str(num_nnis)+' nnis', end='', file=sys.stderr)
        #print(PL)
        
    print(' done', file=sys.stderr)
    print(tree)
    print(PL)
    return tree,PL_new

if __name__ == '__main__':
    import argparse
    
    GTYPE3 = np.array(('RR','RA','AA'))
    GTYPE10 = np.array(('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'))
    DELTA = 1e-4

    parser = argparse.ArgumentParser()
    subp = parser.add_subparsers(metavar='<command>', help='sub-commands')

    #nbjoin uses neighbor_main, read_vcf, make_base_prior (normalize_PL), make_mut_matrix (phred2p, gtype_distance), make_D (pairwise_diff, normalize2d_PL, phred2p), init_star_tree, neighbor_joining
    #init_tree, populate_tree_PL, calc_mut_likelihoods (p2phred), recursive_NNI (nearest_neighbor_interchange, update_PL, score), recursive_reroot (reroot, update_PL)
    parser_nbjoin = subp.add_parser('nbjoin', help='neighbor-joining')
    parser_nbjoin.add_argument('vcf', metavar='<vcf>', type=str, help='input vcf/vcf.gz file, "-" for stdin')
    parser_nbjoin.add_argument('output', metavar='output', type=str, help='output basename')
    parser_nbjoin.add_argument('-m', metavar='INT', dest='mu', type=int, default=80, help='mutation rate in Phred scale, default 80')
    parser_nbjoin.add_argument('-e', metavar='INT', dest='het', type=int, default=30, help='heterozygous rate in Phred scale, default 30')
    parser_nbjoin.add_argument('-v', metavar='INT', dest='min_ev', type=int, default=60, help='minimum evidence in Phred scale for a site to be considered, default 60')
    parser_nbjoin.set_defaults(func=neighbor_main)
    
    try:
        args = parser.parse_args()
        args.func(args)
    except KeyboardInterrupt:
        sys.exit(1)