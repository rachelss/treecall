#!/usr/bin/env bash
#script; num samples; cov; ref; folder containing simulations

mu=$1
if [ "$mu" = "" ]
then
    echo "enter mu"
    read mu
fi

#move previous tree scores
scorefiles=*/x*/x*.treecall.scores.txt
scorefiles=($scorefiles)
if [ "${#scorefiles[@]}" -gt 0 ]; then
    echo "Enter the previous mu"
    read old_mu
fi
if [ "${#scorefiles[@]}" -gt 0 ]; then
    for f in "${scorefiles[@]}"; do
        mv $f $( echo ${f} | sed s/\.treecall\.scores.txt/mu${old_mu}\.treecall\.scores\.txt/ )
    done
fi

bash sim_change_mu.sh 5 5 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 7 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 10 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 15 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 20 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 30 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 40 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 5 50 ref/chr22_20-21M.fa ms5i100s35a $mu
bash sim_change_mu.sh 10 5 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 7 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 10 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 15 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 20 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 30 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 40 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 10 50 ref/chr22_20-21M.fa ms10i100s35a $mu
bash sim_change_mu.sh 20 5 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 7 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 10 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 15 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 20 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 30 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 40 ref/chr22_18-23M.fa ms20i100s35a $mu
bash sim_change_mu.sh 20 50 ref/chr22_18-23M.fa ms20i100s35a $mu