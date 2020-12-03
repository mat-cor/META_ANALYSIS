#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value"}
NR >1 {print $a["chromosome"],$a["base_pair_location"],$a["other_allele"],$a["effect_allele"],0.5,2,$a["beta"],$a["standard_error"],$a["p_value"]}' | \
bgzip -@4 > $1.formatted.gz
