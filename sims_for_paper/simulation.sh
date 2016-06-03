#!/usr/bin/env bash
#bcftools 1.3; Using htslib 1.3

N=$1
cov=$2
ref=$3

if [ "$ref" = "" ]
then
    echo "usage: $0 <N> <cov> <ref>"
    exit 0
fi

samples=`seq 0 $N`
dir="x$cov"
treecall = ../treecall.py
pyfilter = find_polymorphic_sites.py
mkphylip = vcf2seq.py
#phylip tools must be installed and in path dnaml, dnacomp, dnadist, neighbor

which dnaml &>/dev/null
    [ $? -eq 0 ] || { echo "Phylip must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }
which samtools &>/dev/null
    [ $? -eq 0 ] || { echo "Samtools must be installed to run. The installation folder must be in your path. Aborting."; exit 1; }

#commented out b/c using previous sims
## make BAMs
#for s in $samples
#do
#    mkdir -p $dir/$s && dwgsim -c 0 -e 0.005-0.01 -E 0.005-0.01 -1 100 -2 100 -d 350 -s 30 -C $cov -r 0 -m var/$s.variants.txt $ref $dir/$s/$s
#    rm -f $dir/$s/$s.bfast.*
#    gzip -f $dir/$s/$s.bwa.read1.fastq
#    gzip -f $dir/$s/$s.bwa.read2.fastq
#
#    header="@RG\\tID:$s\\tSM:s$s"
#    bwa mem -I350,35 -R$header $ref $dir/$s/$s.bwa.read1.fastq.gz $dir/$s/$s.bwa.read2.fastq.gz | samtools1 sort -o $dir/$s/$s.bam -T $dir/$s/$s -
#    samtools1 index $dir/$s/$s.bam
#done
#echo "BAMs ready"
#
find $dir -name '[1-9]*.bam' | sort-alt -N > $dir/bam.list
find $dir -name '[0-9]*.bam' | sort-alt -N > $dir/bam0.list

# make BCFs
samtools mpileup -f $ref -b $dir/bam.list -C50 -q 13 -Q 13 --ff 0x704 -t AD,ADF,ADR -g -o $dir/$dir.bcf  
samtools mpileup -f $ref -C50 -q 13 -Q 13 --ff 0x704 -t AD,ADF,ADR -g -o $dir/0/0.bcf $dir/0/0.bam
echo "BCFs ready"

# bcfcall
bcftools call -cv -o $dir/$dir.snp.vcf.gz -O z $dir/$dir.bcf  #consensus call, variants only, output compressed vcf
bcftools call -cv -o $dir/0/0.snp.vcf.gz -O z $dir/0/0.bcf
echo "bcfcall done"

# prepare for treecall - output of bcftools view not compatible with vcf.Reader in pyfilter
bcftools view -v snps -i 'ADF[0]>0 & ADF[1]>0 & ADR[0]>0 & ADR[1]>0' $dir/$dir.bcf | #get just snps; reads on both strands to ensure no weird indels - coudl actually test whether prob or quality is same forward or reverse
python2 $pyfilter - 'AD:2;PL:60' | #pipe to filter to check both alleles have more than min reads AND both het and homo present AND likely more than one geno across all samples
bgzip -c > $dir/$dir.candidate1.vcf.gz  #pipe to zip; Skip lines started with character CHAR. [#]

#fix below to allow for gz
python2 $pyfilter $dir/$dir.candidate1.vcf.gz 'AD4:1' | bgzip -c > $dir/$dir.candidate2.vcf.gz #filter for at least F or R has 1 for both refs and alts
#not sure this is dif from filtering just AD

# infer tree
#fix below to allow for gz
#python2 $treecall part -m 60 -v 60 -e 30 $dir/$dir.candidate2.vcf.gz $dir/$dir.treecall-pt
python2 $treecall nbjoin -m 60 -v 60 -e 30 $dir/$dir.candidate2.vcf.gz $dir/$dir.treecall-nj

# treecall genotype
bcftools view $dir/$dir.bcf | python2 $pyfilter - 'AD:2' | python2 $treecall gtype -t $dir/$dir.treecall-pt.nwk -m 60 -e 30 - $dir/$dir.tc.txt
awk '$5>0.5' $dir/$dir.tc.txt > $dir/$dir.tc.p50.txt

bcftools view $dir/$dir.bcf | python2 $pyfilter - 'AD:2' | python2 $treecall gtype -t ms.nwk -m 60 -e 30 - $dir/$dir.ms.txt
awk '$5>0.5' $dir/$dir.ms.txt > $dir/$dir.ms.p50.txt
echo "treecall done"

# phylip
mkdir -p phylip
bcftools view -v snps $dir/$dir.snp.vcf.gz | python2 $mkphylip - phylip > phylip/$dir.snp.phylip
cd phylip
rm -f outtree outfile
echo "$dir.snp.phylip" > phylip_inputs.list
echo "Y" >> phylip_inputs.list
echo "$dir.snp.dist" > dist_inputs.list
echo "Y" >> dist_inputs.list
# --- dnaml --- #
cat phylip_inputs.list | dnaml
paste -s outtree | perl -lne 's/s(\d+)/$1-1/ge; s/:[0-9.-]+/:1/g; s/\s//g; s/\[[^\[\]]*\]//g; s/;/;\n/g; print;' > $dir.ml.tree
head -1 $dir.ml.tree > $dir.dnaml.nwk
rm -f outtree outfile
# --- dnacomp --- #
cat phylip_inputs.list | dnacomp
paste -s outtree | perl -lne 's/s(\d+)/$1-1/ge; s/:[0-9.-]+/:1/g; s/\s//g; s/\[[^\[\]]*\]//g; s/;/;\n/g; print;' > $dir.comp.tree
head -1 $dir.comp.tree > $dir.dnacomp.nwk
rm -f outtree outfile
# --- neighbor --- #
cat phylip_inputs.list | dnadist
mv outfile $dir.snp.dist
cat dist_inputs.list | neighbor
paste -s outtree | perl -lne 's/s(\d+)/$1-1/ge; s/:[0-9.-]+/:1/g; s/\s//g; s/\[[^\[\]]*\]//g; s/;/;\n/g; print;' > $dir.nj.tree
head -1 $dir.nj.tree > $dir.neighbor.nwk
rm -f outtree outfile
cd ..
echo "phylip done"
