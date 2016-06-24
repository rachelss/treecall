import sys
sys.path.append('../')
from tree_est import *
from ete2 import Tree
import numpy as np

def test_read_vcf():
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)
    
    assert len(variants) == 65, "Finding "+len(variants)+" instead of 65 expected"
    
    expectedvarpos = ['1504', '2229', '3689', '15129', '22495', '25009', '28659', '37734', '40100', '41673',
                  '42316', '42349', '45427', '48245', '49676', '64953', '135736', '137570', '143975', '154547',
                  '156846', '183819', '245368', '254022', '271559', '272320', '294334', '299078', '301270',
                  '301574', '338665', '338956', '341437', '344815', '348072', '349640', '350082', '353617',
                  '354515', '358983', '369469', '380869', '382963', '385502', '385929', '387317', '389904',
                  '401804', '402190', '402945', '408259', '408531', '412783', '421695', '422262', '445318',
                  '460513', '468460', '474204', '475512', '478325', '488270', '495948', '499773', '508919']
    
    assert expectedvarpos == list(variants[:,1]), 'Not getting the correct correct variant positions'
    
    ads0 = np.array([[ 4,  0],     [ 5,  7],     [13,  0],     [ 2,  0],     [ 9,  0],     [ 8,  0],     [16,  0],     [ 4,  0],     [13,  0],     [15,  0]])
    assert np.array_equal(ADs[0], ads0), 'The first set of allelic depths has not been read correctly'
    
    pls0 = np.array([[  0,  12,  68], [ 83,   0,  51], [  0,  39, 177], [  0,   6,  43], [  0,  27, 138], [  0,  24, 136], [  0,  48, 187], [  0,  12,  74], [  0,  39, 178], [  0,  45, 166]])
    assert np.array_equal(PLs[0], pls0)
    
    n_site,n_smpl,n_gtype = PLs.shape
    assert n_site == 65
    assert n_smpl == 10
    assert n_gtype == 3
    
def test_init_tree():
    
    tree = Tree("(1:1,(2:1,(5:1,4:1):0.5):0.5);" )
    
    tree = init_tree(tree)
    
    nid = [node.nid for node in tree.traverse(strategy='postorder')]
    sid = [node.sid for node in tree.traverse(strategy='postorder')]
    
    assert nid == [0, 1, 2, 3, 4, 5, 6]
    assert sid == [[1], [2], [5], [4], [4, 5], [2, 4, 5], [1, 2, 4, 5]]
    
#def test_make_base_prior():
#    GTYPE3 = np.array(('RR','RA','AA'))
#    bp = make_base_prior(30, GTYPE3)
#    exp_bp = np.array([3.0124709, 33.012471, 3.0124709])
#    
#    assert(np.isclose(bp,exp_bp).all()==True), 'Base prior incorrect'

#def test_make_mut_matrix():
#    GTYPE3 = np.array(('RR','RA','AA'))
#    mm,mm0,mm1 = make_mut_matrix(80, GTYPE3)
#    
#    exp_mm = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#    exp_mm0 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
#    exp_mm1 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
#    
#    assert(np.isclose(exp_mm, mm).all()==True)
#    assert(np.isclose(exp_mm0, mm0).all()==True)
#    assert(np.isclose(exp_mm1, mm1).all()==True)

def test_make_D():
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)
    PLs = PLs.astype(np.longdouble)
    D = make_D(PLs)
    
    exp_D0 = np.array([ 0.0, 33.871109, 14.227026, 13.116167, 9.5826616, 35.841949, 11.934257, 11.326775, 13.333431, 43.270412, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    assert(np.isclose(D[0], exp_D0).all()==True)

def test_init_star_tree():
    tree = init_star_tree(10)
    assert str(tree.write()) == '(0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1);'

def test_neighbor_joining():
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)
    PLs = PLs.astype(np.longdouble)
    Da = make_D(PLs)
    n_site,n_smpl,n_gtype = PLs.shape
    tree = init_star_tree(n_smpl)
    internals = np.arange(n_smpl)
    D,tree = neighbor_joining(Da.copy(), tree, internals)
    
    exp_D0 = np.array([ 0.0, 33.871108, 14.227027, 13.116167, 9.5826634, 35.841953, 11.934258, 11.326775, 13.333432, 43.27041, 60.619184, 77.172342, 12.930811, 15.712338, -15.291832,  0.0 , 0.0, 0.0])
    assert(np.isclose(D[0], exp_D0).all()==True), 'Error updating difference matrix'
    
    exp_tree = '((4:-50,(2:-26,(1:-4,(5:5,9:12)1:21)1:103)1:62)1:1,(8:-107,(0:-15,(7:-5,(3:6,6:5)1:14)1:31)1:125)1:1);'
    assert str(tree.write()) == exp_tree, str(tree.write())

def test_init_tree():
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)
    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)  # pairwise differences between samples based only on PLs (should include mutation, but also shouldn't matter)
    tree = init_star_tree(n_smpl)
    internals = np.arange(n_smpl)
    D,tree = neighbor_joining(D.copy(), tree, internals) #haven't checked this; make nj tree and update D given internal nodes; pass copy

    tree = init_tree(tree)
    nids,sids=[],[]
    for node in tree.traverse(strategy='postorder'):
        nids.append(node.nid)
        sids.append(node.sid)
    assert nids == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    assert sids == [[4], [2], [1], [5], [9], [5, 9], [1, 5, 9], [1, 2, 5, 9], [1, 2, 4, 5, 9], [8], [0], [7], [3], [6], [3, 6], [3, 6, 7], [0, 3, 6, 7], [0, 3, 6, 7, 8], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    
def test_populate_tree_PL():    
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)

    GTYPE3 = np.array(('RR','RA','AA'))
    mm,mm0,mm1 = make_mut_matrix_gtype3(80) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)  # pairwise differences between samples based only on PLs (should include mutation, but also shouldn't matter)
    tree = init_star_tree(n_smpl)
    internals = np.arange(n_smpl)
    D,tree = neighbor_joining(D.copy(), tree, internals) #haven't checked this; make nj tree and update D given internal nodes; pass copy

    tree = init_tree(tree)  #tree has nid's (node id) and sid's (list of tip names - sorted)
    tree = populate_tree_PL(tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
    
    PL_var_0_tips = []
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf():
            PL_var_0_tips.append((node.sid,node.PL0[0]))  #get PLs for var0 for each tip
    
    PL_var_0_tips.sort(key=lambda tup: tup[0])  #sort by sample numbers
    PL_var_0_tips = np.array(PL_var_0_tips)[:,1]
#    print(PL_var_0_tips)
#    print(PLs[0,:])
    for i,j in enumerate(PL_var_0_tips):
        assert np.array_equal(PL_var_0_tips[i], PLs[0,:][i])  #compare original PLs to what's in tree
        
def test_calc_mut_likelihoods():    
    input_vcf = 'test_tree.vcf'
    vcffile, variants, ADs, PLs = read_vcf(input_vcf, 60)

    GTYPE3 = np.array(('RR','RA','AA'))
    mm,mm0,mm1 = make_mut_matrix_gtype3(80) # substitution rate matrix, with non-diagonal set to 0, with diagonal set to 0

    PLs = PLs.astype(np.longdouble)
    n_site,n_smpl,n_gtype = PLs.shape

    D = make_D(PLs)  # pairwise differences between samples based only on PLs (should include mutation, but also shouldn't matter)
    tree = init_star_tree(n_smpl)
    internals = np.arange(n_smpl)
    D,tree = neighbor_joining(D.copy(), tree, internals) #haven't checked this; make nj tree and update D given internal nodes; pass copy

    tree = init_tree(tree)  #tree has nid's (node id) and sid's (list of tip names - sorted)
    tree = populate_tree_PL(tree, PLs, mm0, 'PL0')  #tree has PLs for no mutation at tips and nodes
    tree = calc_mut_likelihoods(tree, mm0, mm1)  #attach PLm to each node (not tips!)
    
    