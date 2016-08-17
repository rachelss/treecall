import sys
import string
import random
from Bio import SeqIO
import os
from collections import Counter

def get_pos(poslist):
    """list of positions with mutations"""
    poslist = poslist[1:]
    poslist = [float(p) for p in poslist]
    #change duplicate mutations
    p_counts = Counter(poslist)
    for k,v in p_counts.iteritems():
        while v > 1:
            v = v-1
            pos_index = poslist.index(k)  #get position that is duplicate
            poslist[pos_index+v] = k + (v*0.00001)     #change last of duplicate positions by adding a little        
    poslist = [int(p*1000000) for p in poslist]
    assert len(poslist) > len(set(poslist)), 'duplicate mutations at a position'
    return poslist

def get_refbases(seq,poslist):
    """list of reference bases for each variant position"""
    refbases = [seq[p] for p in poslist]
    return refbases

def get_altbases(refbases):
    """ list of IUPAC codes for hets"""
    altbases=[]
    bases = ['A','G','C','T']
    IUPAC = { ('A','C'):'M', ('A','G'):'R', ('A','T'):'W', ('C','G'):'S', ('C','T'):'Y', ('G','T'):'K' }

    for r in refbases:
        b=bases[:]
        b.remove(r)
        newbase = random.choice(b)
        baselist = tuple(sorted([r,newbase]))
        altbases.append(IUPAC[baselist])
    return altbases

def get_strands(num_sites):
    strands = [random.randint(1,2) for i in range(num_sites)]
    return strands

###############
folder = os.path.dirname(sys.argv[1])

msfile = open(sys.argv[1],'r')
msoutput = msfile.readlines()
cmd = msoutput[0].split()
num_samp = int(cmd[1])
num_sites = int(cmd[4])

pos = get_pos(msoutput[6].split())

handle = open('ref/chr22_20-21M.fa', "rU")
ref = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
handle.close()

refbases = get_refbases(ref['chr22'].seq,pos)
refbases = [r.upper() for r in refbases]
altbases = get_altbases(refbases)
strands = get_strands(num_sites)

for i in range(num_samp):
    v = str(i+1)
    dwgsimfile = open(folder+'/var/'+v+'.variants.txt','w')
    linenum = i-num_samp
    line = msoutput[linenum]
    for j in range(num_sites):
        if int(line[j]) == 1:
            #chr    #pos    #refbase    #altbase/iupac  #strand
            dwgsimfile.write("chr22\t"+str(pos[j]+1)+"\t"+refbases[j]+"\t"+altbases[j]+"\t"+str(strands[j])+"\n")
    dwgsimfile.close()         

treefile = open(sys.argv[1].replace('.output','.nwk'),'w')
treefile.write(msoutput[4])
treefile.close()