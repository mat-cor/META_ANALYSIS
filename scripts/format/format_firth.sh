#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value"}
NR >1 {print $a["chromosome"],$a["base_pair_location"],$a["other_allele"],$a["effect_allele"],0.5,2,$a["beta"],$a["standard_error"],$a["p_value"]}' | \
bgzip -@4 > $1.formatted.gz


#CHROMPOSIDREFALTA1FIRTH?TESTOBS_CTORLOG(OR)_SEZ_STATPERRCODE
1634553rs368347679GAANADD2880.6292880.678478-0.6826550.494825.
1861945rs111739932TCCNADD2911.89960.6976040.9197840.357686.
1866281rs12132974CTTNADD2911.899510.6975980.9197220.357718.
