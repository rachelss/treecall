#!/usr/bin/env bash
#bcftools 1.3; Using htslib 1.3

N=$1  #number of samples
cov=$2  #coverage
ref=$3  #reference for simulation and alignment
basefolder=$4  #folder containing simulations

if [ "$basefolder" = "" ]
then
    echo "usage: $0 <N> <cov> <ref> <basefolder>"
    exit 0
fi

samples=`seq 0 $N`
dir="${basefolder}/x${cov}"
treecall=../treecall.py
pyfilter=find_polymorphic_sites.py
mkphylip=vcf2seq.py

#phylip must be installed and in path for dnacomp
which dnacomp &>/dev/null
    [ $? -eq 0 ] || { echo "Phylip must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }
which raxml &>/dev/null
    [ $? -eq 0 ] || { echo "Raxml must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }
which samtools &>/dev/null
    [ $? -eq 0 ] || { echo "Samtools must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }
which dwgsim &>/dev/null
    [ $? -eq 0 ] || { echo "dwgsim must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }
which bwa &>/dev/null
    [ $? -eq 0 ] || { echo "bwa must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }

# make BAMs
for s in $samples; do
    if [[ $s -eq 0 ]]; then
        mkdir -p $dir/$s && dwgsim -c 0 -e 0.005-0.01 -E 0.005-0.01 -1 100 -2 100 -d 350 -s 30 -C $cov -r 0 $ref $dir/$s/$s
    else
        mkdir -p $dir/$s && dwgsim -c 0 -e 0.005-0.01 -E 0.005-0.01 -1 100 -2 100 -d 350 -s 30 -C $cov -r 0 -m ${basefolder}/var/$s.variants.txt $ref $dir/$s/$s
    fi
    rm -f $dir/$s/$s.bfast.*
    gzip -f $dir/$s/$s.bwa.read1.fastq
    gzip -f $dir/$s/$s.bwa.read2.fastq

    header="@RG\\tID:$s\\tSM:s$s"
    bwa mem -I350,35 -R$header $ref $dir/$s/$s.bwa.read1.fastq.gz $dir/$s/$s.bwa.read2.fastq.gz | samtools sort -o $dir/$s/$s.bam -T $dir/$s/$s -
    samtools index $dir/$s/$s.bam
done
echo "BAMs ready"

find $dir -name '[1-9]*.bam' | sort -V > $dir/bam.list
find $dir -name '[0-9]*.bam' | sort -V > $dir/bam0.list

# make BCFs
samtools mpileup -f $ref -b $dir/bam.list -C50 -q 13 -Q 13 --ff 0x704 -t AD,ADF,ADR -g -o ${dir}/x${cov}.bcf  #output is PL:ADF:ADR:AD for all sites 
samtools mpileup -f $ref -C50 -q 13 -Q 13 --ff 0x704 -t AD,ADF,ADR -g -o $dir/0/0.bcf $dir/0/0.bam
echo "BCFs ready"

# bcfcall
bcftools call -cv -o $dir/x${cov}.snp.vcf.gz -O z ${dir}/x${cov}.bcf   #consensus call, variants only, output compressed vcf
bcftools call -cv -o $dir/0/0.snp.vcf.gz -O z $dir/0/0.bcf
echo "bcfcall done"

#get just snps; reads on both strands to ensure no weird indels - coudl actually test whether prob or quality is same forward or reverse
#bcftools view -v snps -i 'ADF[0]>0 & ADF[1]>0 & ADR[0]>0 & ADR[1]>0' $dir/x${cov}.bcf | sed 's/Number=[A-Z]/Number=./' | sed 's/,Version=\"3\"//' > $dir/x${cov}.var.vcf
#not filtering by F and R strands due to sim, which puts mutation on a strand not chr
bcftools view -v snps $dir/x${cov}.bcf | sed 's/Number=[A-Z]/Number=./' | sed 's/,Version=\"3\"//' > $dir/x${cov}.var.vcf
# prepare for treecall - output of bcftools view not compatible with vcf.Reader in pyfilter
python2 $pyfilter $dir/x${cov}.var.vcf 'AD:2;PL:60'
# > $dir/x${cov}.candidate1.vcf #pipe to filter to check both alleles have more than min reads AND both het and homo present AND likely more than one geno across all samples
#bgzip -c > $dir/x${cov}.candidate1.vcf.gz  #pipe to zip; Skip lines started with character CHAR. [#]

#fix below to allow for gz
#python2 $pyfilter $dir/x${cov}.var.vcf.vcf 'AD4:1'
# | bgzip -c > $dir/x${cov}.candidate2.vcf.gz #filter for at least F or R has 1 for both refs and alts
#not sure this is dif from filtering just AD

# --- treecall infer tree --- #
python2 $treecall nbjoin -m 60 -v 60 -e 30 $dir/x${cov}.var.vcf.vcf $dir/x${cov}.treecall
besttree=$(sort -n -k 2 $dir/x${cov}.treecall.scores.txt | head -1 | cut -f 1 -d ' ')
sed 's/s//g' <$dir/x${cov}.treecall.${besttree}names.tre >$dir/x${cov}.treecall_num.tre
echo "treecall done"

# --- treecall genotyping on best tree--- #
bcftools view $dir/x${cov}.bcf | sed 's/Number=[A-Z]/Number=./' | sed 's/,Version=\"3\"//' > $dir/x${cov}.vcf
python2 $pyfilter $dir/x${cov}.vcf 'AD:2'
python2 $treecall gtype -t $dir/x${cov}.treecall.${besttree}names.tre -m 60 -e 30 $dir/x${cov}.vcf.vcf $dir/x${cov}.tc.txt
awk '$5>0.5' $dir/x${cov}.tc.txt > $dir/x${cov}.tc.p50.txt
echo "genotypes estimated on treecall tree"

# --- treecall genotyping on ms tree--- #
sed 's/\([0-9][0-9]*:\)/s\1/g' <$(dirname $dir)/ms.nwk > $(dirname $dir)/ms2.nwk  #match tree names and vcf sample names by adding s in front of numbers
python2 $treecall gtype -t $(dirname $dir)/ms2.nwk -m 60 -e 30 $dir/x${cov}.vcf.vcf $dir/x${cov}.ms.txt
awk '$5>0.5' $dir/x${cov}.ms.txt > $dir/x${cov}.ms.p50.txt
echo "genotypes estimated on ms tree"

# --- get fixed genotypes - var sites --- #
mkdir -p phylip
bcftools view -v snps $dir/x${cov}.snp.vcf.gz > $dir/x${cov}.snp.vcf
python2 $mkphylip $dir/x${cov}.snp.vcf

# -- dnacomp -- #
echo "$dir/x${cov}.snp.vcf.alignment.phylip" > phylip/phylip_inputs${basefolder}${cov}.list
echo "Y" >> phylip/phylip_inputs${basefolder}${cov}.list

# --- raxml --- #
#ascertainment bias correction; ASC_GTRCAT -V = plain GTR
raxml -s "$dir/x${cov}.snp.vcf.alignment.phylip" -n out${basefolder}${cov} -m ASC_GTRCAT -V --asc-corr=lewis -T 4 -p $RANDOM
sed 's/s//g' <RAxML_bestTree.out${basefolder}${cov} >$dir/x${cov}.ml_num.tre
rm RAxML*out${basefolder}${cov}*


# --- neighbor --- #
#echo "$dir.snp.dist" > dist_inputs.list
#echo "Y" >> dist_inputs.list
#cat phylip_inputs.list | dnadist
#mv outfile $dir.snp.dist
#cat dist_inputs.list | neighbor
#paste -s outtree | perl -lne 's/s(\d+)/$1-1/ge; s/:[0-9.-]+/:1/g; s/\s//g; s/\[[^\[\]]*\]//g; s/;/;\n/g; print;' > $dir/x${cov}.nj.tree
#head -1 $dir/x${cov}.nj.tree > $dir/x${cov}.neighbor.tre
#rm -f outtree outfile

#don't use dnaml - no way to correct for ascertainment bias
# --- dnaml --- #
#cat phylip/phylip_inputs.list | dnaml
#paste -s outtree | perl -lne 's/s(\d+)/$1-1/ge; s/:[0-9.-]+/:1/g; s/\s//g; s/\[[^\[\]]*\]//g; s/;/;\n/g; print;' > $dir/x${cov}.ml.tree
#head -1 $dir/x${cov}.ml.tree > $dir/x${cov}.ml.tre
#rm -f outtree outfile
#echo "phylip done"

#-N Number of same-patient samples to analyse (required)
n=`cat $sim/x$cov/bam0.list | wc -l`
#--fasta Reference fasta file
#--medianT Median depth of coverage in tumour samples
#--medianN Median depth of coverage in normal sample
#--bam List of bam files. Place normal bam file FIRST!
#-d Minimum required depth (in normal and as an average across all samples)
d=`echo $cov/2 | bc`  #divides cov / 2
if [ $d -gt 5 ]; then 
    d=5
fi
#-f name of VCF output file
multiSNV -N$n --fasta ref/chr22_20-21M.fa --medianN $cov --medianT $cov -f $sim/x$cov/x$cov.multiSNV.vcf -d $d --bam $(cat $sim/x$cov/bam0.list |tr '\n' ' ')


#----remove bam and vcf----#
rm ${dir}/*/*bam
#rm ${dir}/*vcf