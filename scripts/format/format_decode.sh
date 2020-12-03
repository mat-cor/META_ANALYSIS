#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {FS=" "; OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value"}
NR >1 {print $a["Chrom"],$a["Pos_b38"],$a["Other_allele"],$a["Effect_allele"],$a["EA_freq_PC"]/100,$a["Info"],$a["Beta"],$a["Stderr"],$a["Pval"]}' | \
bgzip -@4 > $1.formatted.gz
