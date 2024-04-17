#!/usr/bin/env bash


set -eu


#######################################
#input genotype
spresto_input=../2_genotype/s969
gustomom_input=../2_genotype/m1079_snpQC

#input data

dir=QC
spresto=$dir/spresto
gusto=$dir/gusto
merged2=$dir/merged2
merged=$dir/merged
mergedQC=$dir/mergedQC
mergedQC_pruned=$dir/mergedQC_pruned

#module load plink/2.00a2.3
if [ ! -e $spresto.bim ]; then
plink2 	--pfile $spresto_input \
	--keep-allele-order \
	--rm-dup  \
	--geno 0.01 \
	--maf 0.05 \
	--hwe 1e-6 \
	--snps-only just-acgt \
	--extract SNPlist.txt \
	--exclude merged-merge.missnp \
	--keep spresto934.txt \
	-output-chr 26 \
	--make-bed --out $spresto
fi

if [ ! -e $gusto.bim ] ; then
plink2 	--pfile $gustomom_input \
	--keep-allele-order \
 	--rm-dup  \
	--geno 0.01 \
	--maf 0.05 \
	--hwe 1e-6 \
	--snps-only just-acgt \
	--extract SNPlist.txt \
	--exclude merged-merge.missnp \
	--keep gusto434.txt \
	-output-chr 26 \
	--make-bed --out $gusto
fi

x=$(cat SNPlist.txt |wc -l)
y=$(cat $spresto.bim |wc -l)
z=$(cat $gusto.bim |wc -l)

echo -e "\t  MetS SNPs   		= $x "
echo -e "\t  MetS SNPs in spresto	= $y "
echo -e "\t  MetS SNPs in gusto   	= $z "

#--exclude $dir/merged-merge.missnp \
#$dir/merged-merge.missnp will automatically
#generated when fail in the first run
#module load plink/1.9
if [ ! -e $merged.bim ] ; then
grep -wf $spresto.bim $gusto.bim > $dir/common.snplist 
plink 	--bfile $spresto  \
	--bmerge $gusto   \
	--keep-allele-order \
	--exclude merged-merge.missnp \
	--extract $dir/common.snplist \
	--make-bed --out $merged
fi


x=$(cat $merged.bim |wc -l)
y=$(cat $merged.fam |wc -l)
echo -e "\t  MetS SNPs merged   	= $x "
echo -e "\t  MetS samples merged (934+434+8822)    	= $y "


###1. extract 
###spresto934=spresto934.txt
###MetSSNPs=SNPlist.txt
###2. remove MHC
###MHC=$dir/MHC.snplist
#module load plink/1.9
if [ ! -e $dir/MHC.snplist ] ; then
plink   --bfile $spresto \
        --chr 6 \
        --from-bp 25477797 \
        --to-bp 36448354   \
        --write-snplist    \
        --out $dir/MHC
fi

x=$(cat $dir/MHC.snplist |wc -l)
echo -e "\t  SNPs in MHC   = $x "

if [ ! -e $mergedQC.bim ]; then
plink \
	--bfile $merged \
	--keep-allele-order \
	--geno 0.01 \
	--maf 0.05 \
	--hwe 1e-10 \
	--snps-only just-acgt \
	--exclude $dir/MHC.snplist \
	--make-bed \
	--out $mergedQC 
fi
x=$(cat $mergedQC.bim |wc -l)
y=$(cat $mergedQC.fam |wc -l)
echo -e "\t  MetS SNPs in merged after QC and removal of MHC   = $x "
echo -e "\t  MetS samples in merged after QC and removal of MHC   = $y "


if [ ! -e $mergedQC_pruned.prune.in ]; then
plink \
        --bfile $mergedQC \
        --indep-pairwise 200 100 0.1 \
        --out $mergedQC_pruned
fi

y=$(cat $mergedQC_pruned.prune.in|wc -l)
z=$(cat $mergedQC_pruned.prune.out|wc -l)
echo -e "\t  prune in = $y"
echo -e "\t  prune out = $z"

if [ ! -e $mergedQC_pruned.bim ]; then
plink   --bfile $mergedQC \
        --extract $mergedQC_pruned.prune.in  \
        --make-bed \
        --out $mergedQC_pruned
fi

x=$(cat $mergedQC_pruned.bim |wc -l)
y=$(cat $mergedQC_pruned.fam |wc -l)
echo -e "\t  final SNPs =$x "
echo -e "\t  final Samples =$y "

if [ ! -e $mergedQC_pruned.traw ]; then 
echo "extract the dosage and transpose it"
plink 	--bfile $mergedQC_pruned  \
	--export A-transpose \
	--out $mergedQC_pruned
fi

echo -e "done"


#Rscript to get the dosage
#####################
rm(list=ls())
library(dplyr)
library(data.table)


dosagefile="QC/mergedQC_pruned.traw"
output="s934_g434_snplist_data.txt"
output_transpose="s934_g434_snplist_data_transposed.txt"

mydata1<- fread(dosagefile, header=T, sep="\t", check.names=F)
colnames(mydata1)[1:10]
colnames(mydata1) <- gsub("_.*", "", colnames(mydata1)); 
mydata1$SNP_ALT <- paste0(mydata1$SNP, "_", mydata1$COUNTED)
dim_y =dim(mydata1)[2]

mycol<- c("SNP_ALT", colnames(mydata1)[7:dim_y])
mydata1<- mydata1 %>% select(all_of(mycol))
dim(mydata1)
mydata1[1:3, 1:3]

#transpose data frame
mydata_t <- transpose(mydata1)
#redefine row and column names
rownames(mydata_t) <- colnames(mydata1)
colnames(mydata_t) <- rownames(mydata1)
fwrite(mydata_t, output_transpose, col.names=F,row.names=T, sep="\t", quote=F)





