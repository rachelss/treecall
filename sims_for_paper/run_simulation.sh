#!/usr/bin/env bash

rm treecomp.txt

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
                sed 's/s//g' <outtree >$dir/x${cov}.dnacomp_num.tre
                rm -f outfile
                rm -f outtree
                
                for testtree in $dir/x${cov}.ml_num.tre $dir/x${cov}.dnacomp_num.tre $dir/x${cov}.treecall_num.tre; do
                    python2 ../treecall.py compare -t $testtree -r ${basefolder}/ms.nwk >> treecomp.txt
                done
            done
        done
    done
done
echo "dnacomp done"