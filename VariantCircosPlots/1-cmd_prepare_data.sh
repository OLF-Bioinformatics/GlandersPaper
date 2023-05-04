#!/bin/bash

# for this commend you need to have vcftools

DIR_DATA=../data
DIR_DATA_CIR=../files

if [ ! -d $DIR_DATA_CIR ]; then
	mkdir -v $DIR_DATA_CIR
fi
 
#====================== SNPs ======================
VCF_SNP_NO_HET=$DIR_DATA/merged.snps.vcf
vcftools --vcf $VCF_SNP_NO_HET --SNPdensity 20000 --out $DIR_DATA_CIR/SNP_density_NO_HET 

SNP_IN=$DIR_DATA_CIR/SNP_density_NO_HET.snpden
SNP_OUT=$DIR_DATA_CIR/snp_density_no_het.txt
echo "snp file: $SNP_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $SNP_IN > $SNP_OUT

# strain specific SNPs
STRAINS="10229 11 2000031063 2002721276 2002734299 2002734306 23344 34 6 ATCC23344_CFIA Bahrain1 BMQ Bogor FDAARGOS_585 FDAARGOS_588 FDAARGOS_589 FDAARGOS_590 FMH23344 FMH India86-567-2 JHU KC_1092 Muk NCTC_10247_1 NCTC_10247 SAVP1 Turkey10 Turkey1 Turkey2 Turkey3 Turkey4 Turkey5 Turkey6 Turkey7 Turkey8 Turkey9 Zagreb"
for strain in $STRAINS;
do
	echo $strain
	strain_vcf=$DIR_DATA_CIR/$strain\_SNP_no_Het.vcf
	bcftools view -s $strain $VCF_SNP_NO_HET | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 20000 --out $DIR_DATA_CIR/$strain\_SNP_density_NO_HET
	SNP_IN=$DIR_DATA_CIR/$strain\_SNP_density_NO_HET.snpden
	SNP_OUT=$DIR_DATA_CIR/$strain\_snp_density_no_het.txt
	echo "snp file: $SNP_OUT"
	awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $SNP_IN > $SNP_OUT
done

#====================== INDELs ======================
VCF_INDEL_NO_HET=$DIR_DATA/merged.indels.vcf
vcftools --vcf $VCF_INDEL_NO_HET --SNPdensity 20000 --out $DIR_DATA_CIR/INDEL_density_NO_HET

INDEL_IN=$DIR_DATA_CIR/INDEL_density_NO_HET.snpden
INDEL_OUT=$DIR_DATA_CIR/indel_density_no_het.txt
echo "INDEL file: $INDEL_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $INDEL_IN > $INDEL_OUT

for strain in $STRAINS;
do
	echo $strain
	strain_vcf=$DIR_DATA_CIR/$strain\_INDEL_no_Het.vcf
	bcftools view -s $strain $VCF_INDEL_NO_HET | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 20000 --out $DIR_DATA_CIR/$strain\_INDEL_density_NO_HET
	INDEL_IN=$DIR_DATA_CIR/$strain\_INDEL_density_NO_HET.snpden
	INDEL_OUT=$DIR_DATA_CIR/$strain\_indel_density_no_het.txt
	echo "INDEL file: $INDEL_OUT"
	awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $INDEL_IN > $INDEL_OUT
done

#====================== SVs ======================
VCF_SV_IN=$DIR_DATA/merged.sv.vcf
vcftools --vcf $VCF_SV_IN --SNPdensity 20000 --out $DIR_DATA_CIR/SV_density

SV_IN=$DIR_DATA_CIR/SV_density.snpden
SV_OUT=$DIR_DATA_CIR/sv_density.txt
echo "sv file: $SV_OUT"
awk 'BEGIN{FS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $SV_IN > $SV_OUT

for strain in $STRAINS;
do
	echo $strain
	strain_vcf=$DIR_DATA_CIR/SV_$strain.vcf
	bcftools view -s $strain $VCF_SV_IN | bcftools view -i 'GT=="1/1"' > $strain_vcf
	vcftools --vcf $strain_vcf --SNPdensity 20000 --out $DIR_DATA_CIR/SV_$strain
	SV_IN=$DIR_DATA_CIR/SV_$strain.snpden
	SV_OUT=$DIR_DATA_CIR/SV_$strain.txt
	echo "SV file: $SV_OUT"
	awk 'BEGIN{OFS="\t"}{if($1!="CHROM"){gsub("chr", "mm", $1); print $1,$2,$2+2e4,$3}}' $SV_IN > $SV_OUT
done
