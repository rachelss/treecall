#!/usr/bin/env bash
#bcftools 1.3; Using htslib 1.3

N=$1  #number of samples
cov=$2  #coverage
ref=$3  #reference for simulation and alignment
basefolder=$4  #folder containing simulations
mu=$5

if [ "$mu" = "" ]
then
    echo "usage: $0 <N> <cov> <ref> <basefolder> <mu>"
    exit 0
fi

dir="${basefolder}/x${cov}"
treecall=../treecall.py

python2 $treecall nbjoin -m $mu -v 60 -e 30 $dir/x${cov}.var.vcf.vcf.vcf $dir/x${cov}mu${mu}.treecall  #tree name has mu in it - don't overwrite original mu60
besttree=$(sort -n -k 2 $dir/x${cov}mu${mu}.treecall.scores.txt | head -1 | cut -f 1 -d ' ')

# treecall genotype
python2 $treecall gtype -t $dir/x${cov}mu${mu}.treecall.${besttree}.tre -m $mu -e 30 $dir/x${cov}.vcf.vcf $dir/x${cov}mu${mu}.tc.txt
awk '$5>0.5' $dir/x${cov}mu${mu}.tc.txt > $dir/x${cov}mu${mu}.tc.p50.txt

mv $dir/x${cov}.treecall.scores.txt $dir/x${cov}mu${mu}.treecall.scores.txt

echo "treecall done"