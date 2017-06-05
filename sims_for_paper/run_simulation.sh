#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo "usage: $0 <processors>"
    exit
fi

bwa index ref/chr22_20-21M.fa

rm sim_list.txt
for num_samp in 5 10 20; do
    for seg_sites in 100 500 1000; do
        for r in {1..10}; do
            #list of sims to run
            #eg bash simulation.sh 5 5 ref/chr22_20-21M.fa ms5i100s35a
            #script; num samples; cov; ref; folder containing simulations
            for cov in 5 7 10 15 20 30 40 50; do
                echo bash simulation.sh $num_samp $cov ref/chr22_20-21M.fa ms${num_samp}i${seg_sites}s${r}r >> sim_list.txt
            done
            
            #run ms
            mkdir -p ms${num_samp}i${seg_sites}s${r}r/var
            ms $num_samp 1 -s $seg_sites -T > ms${num_samp}i${seg_sites}s${r}r/ms.output
            
            #rearrange data into simulation input
            python2 rearrange_ms_dwgsim.py ms${num_samp}i${seg_sites}s${r}r/ms.output

        done
    done
done

cat sim_list.txt | parallel -j $1

rm treecomp.txt
# -- dnacomp and compare trees-- #
for num_samp in 5 10 20; do
    for seg_sites in 100 500 1000; do
        for cov in 5 7 10 15 20 30 40 50; do
            for r in {1..10}; do
                basefolder=ms${num_samp}i${seg_sites}s${r}r
                dir="${basefolder}/x${cov}"
                
                cat phylip/phylip_inputs${basefolder}${cov}.list | dnacomp
                cat outtree | tr '\n' '@' | sed 's/,@/,/g' | tr '@' '\n' | head -1 | sed 's/s//g' | sed 's/\[.*\]//' >$dir/x${cov}.dnacomp_num.tre  #get rid of weird newlines, get first tree, just num names
                rm -f outfile
                rm -f outtree
                
                for testtree in $dir/x${cov}.ml_num.tre $dir/x${cov}.dnacomp_num.tre $dir/x${cov}.treecall_num.tre; do
                    python2 ../treecall.py compare -t $testtree -r ${basefolder}/ms.nwk >> treecomp.txt
                done
            done
        done
    done
done
echo "dnacomp and tree comparisons done"

sed 's/ms//g' <treecomp.txt | sed 's/i/ /g' | sed 's/s/ /g' | sed 's/_num.tre//g' | sed 's/r / /g' | sed 's/\/x[0-9]*\./ /g' | sed 's/r\/x/ /g' | tr '\t' ' ' >treecomp2.txt

# -- plot comparisons of tree success -- #
Rscript plot_treecomp_rf.R 

# -- compare genotypes -- #
rm eval_list.txt
for num_samp in 5 10 20; do
    for seg_sites in 100 500 1000; do
        for r in {1..10}; do
            basefolder=ms${num_samp}i${seg_sites}s${r}r
            
            #get all true variants by combining variants for individuals
            #output looks like
            #chr22   301     T       K
            cat "${basefolder}"/var/*.variants.txt | sort -n -k2 | uniq | cut -f 1,2,3,4 > "${basefolder}"/var/allvars.txt
            
            #get true genotypes
            #output file is lines of
            #chr22   140     G       G       R       G       G       G       G       G       G       G       G
            ./make_true_gt_matrix.py $basefolder $num_samp $seg_sites  #writes to "${basefolder}"/var/true.spgt.txt

            for cov in 5 7 10 15 20 30 40 50; do                               
                echo bash evaluate2.sh $basefolder $num_samp $cov $seg_sites >> eval_list.txt
            done
        done
    done
done
cat eval_list.txt | parallel -j $1
cat */*/gt_comparisons.txt |sort -n|uniq > gt_comparisons2.txt
echo "genotype comparisons done"

# -- plot comparisons of genotyping success -- #
#cat gt_comparisons.txt|sort -n|uniq >gt_comparisons2.txt
Rscript plot_gt.R


#simulations with allelic dropout (ado)
#identical to regular simulations, but
#first add sites that are het for all samples
rm sim_list.txt
for num_samp in 5 10 20; do
    seg_sites=500
    for r in {1..10}; do
        mkdir -p ms${num_samp}i${seg_sites}s${r}rHET/var
        for f in ms${num_samp}i${seg_sites}s${r}r/var/*s.txt; do cp $f ms${num_samp}i${seg_sites}s${r}rHET/var; done  #copy var files
        ./add_hets.py ms${num_samp}i${seg_sites}s${r}rHET/var 50
        ln ms${num_samp}i${seg_sites}s${r}r/ms.nwk ms${num_samp}i${seg_sites}s${r}rHET
        
        for cov in 5 7 10 15 20 30 40 50; do           
            echo bash simulation.sh $num_samp $cov ref/chr22_20-21M.fa ms${num_samp}i${seg_sites}s${r}rHET >> sim_list.txt
        done
    done
done

cat sim_list.txt | parallel -j $1

#now each added site has some probability of looking like a homozygote