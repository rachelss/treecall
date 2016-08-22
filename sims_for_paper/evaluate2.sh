# sample-wise evaluation of different methods

sim=$1
n_sample=$2
cov=$3
nt=$4

if [ "x$cov" = "x" ]
then
    echo "usage: $0 <sim_dir> <n_sample> <cov> <seg sites>"
    exit
fi

>&2 echo "cov=x$cov"

#get root variable sites for filtering results (chr and pos)
bcftools view -v snps $sim/x$cov/0/0.snp.vcf.gz | egrep -v '#' | cut -f1,2 > $sim/x$cov/0/varsites.txt

#---- header ----#
#need to add segsites
echo "nsample   segsites    cov method  precision   recall" > $sim/x$cov/gt_comparisons.txt

#---- compare bcftools est genotypes ----#
t=30
bcftools view -v snps $sim/x$cov/x$cov.snp.vcf.gz | awk -v t=$t '$1~/^#/ || $6>t' > $sim/x$cov/bcftools_est.vcf
./vcf2gt.py $sim/x$cov/bcftools_est.vcf $sim/x$cov/est.spgt.bcf.txt  #infile outfile
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.bcf.txt $nt`  #values for precision and recall
echo "$n_sample $nt $cov    bcftools    $results" >> gt_comparisons.txt

#---- compare bcftools+filter on root ----#
awk 'NR==FNR{a[$2];next} !($2 in a)' $sim/x$cov/0/varsites.txt $sim/x$cov/est.spgt.bcf.txt > $sim/x$cov/est.spgt.bcffilt.txt
#combine $sim/x$cov/0/varsites.txt "${sim}"/var/est.spgt.bcf.txt -f '1,2:-::1,2:-:' -m 4 | sort-alt -k1,1N -k2,2N 
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.bcffilt.txt $nt`
echo "$n_sample $nt $cov    bcftools+filter    $results" >> $sim/x$cov/gt_comparisons.txt

#---- multisnv ----#
cat $sim/x$cov/x$cov.multiSNV.vcf | awk '$1~/^#/ || $7!~/LOW_QUAL/' > $sim/x$cov/x$cov.multiSNV.vcf
./vcf2gt.py $sim/x$cov/x$cov.multiSNV.vcf $sim/x$cov/est.spgt.msnv.txt
cat $sim/x$cov/est.spgt.msnv.txt|cut -f1-3,5- > $sim/x$cov/est.spgt.msnv2.txt  #remove sample 0
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.msnv2.txt $nt`
echo "$n_sample $nt $cov    msnv    $results" >> $sim/x$cov/gt_comparisons.txt


#---- treecall ----#
t=0.5
awk -v t=$t '$5>t' $sim/x$cov/x$cov.tc.p50.txt > $sim/x$cov/treecall_est.txt
./treecall2gt.py $sim/x$cov/treecall_est.txt $sim/x$cov/est.spgt.treecall.txt $n_sample #infile outfile nsamples
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.treecall.txt $nt`
echo "$n_sample $nt $cov    treecall    $results" >> $sim/x$cov/gt_comparisons.txt

#---- treecall+filter ----#
awk 'NR==FNR{a[$2];next} !($2 in a)' $sim/x$cov/0/varsites.txt $sim/x$cov/est.spgt.treecall.txt > $sim/x$cov/est.spgt.treecallfilt.txt
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.treecallfilt.txt $nt`
echo "$n_sample $nt $cov    treecall+filter    $results" >> $sim/x$cov/gt_comparisons.txt

#---- treecall* ----#
awk -v t=$t '$5>t' $sim/x$cov/x$cov.ms.p50.txt > $sim/x$cov/treecall_est_star.txt
./treecall2gt.py $sim/x$cov/treecall_est_star.txt $sim/x$cov/est.spgt.treecallstar.txt $n_sample
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.treecallstar.txt $nt`
echo "$n_sample $nt $cov    treecall*    $results" >> $sim/x$cov/gt_comparisons.txt

#---- treecall*+filter ----#
awk 'NR==FNR{a[$2];next} !($2 in a)' $sim/x$cov/0/varsites.txt $sim/x$cov/est.spgt.treecallstar.txt > $sim/x$cov/est.spgt.treecallstarfilt.txt
#awk -v t=$t '$5>t' $sim/x$cov/x$cov.ms.p50.txt | combine $sim/x$cov/0/varsites.txt - -f '1,2:-::1,2:-:' -m 4 | sort-alt -k1,1N -k2,2N | python $tree2gt - $n_sample > tmp.txt
results=`./compare_gt_matrix.py "${sim}"/var/true.spgt.txt $sim/x$cov/est.spgt.treecallstarfilt.txt $nt`
echo "$n_sample $nt $cov    treecall*filt    $results" >> $sim/x$cov/gt_comparisons.txt
