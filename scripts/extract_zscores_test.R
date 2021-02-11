#!/usr/bin/env Rscript

require(optparse)
require(data.table)
require(R.utils)
require(ggplot2)
require(dplyr)

option_list = list(
  make_option(c("-g", "--gwas_code"), type="character", default=NULL,
              help="gwas code", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL,
              help="chromosome", metavar="character"),
  make_option(c("-r", "--ref_path"), type="character", default=NULL,
              help="ref panel path", metavar="character"),
  make_option(c("-n", "--n_snps"), type="character", default=NULL,
              help="n snps to extract", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

key <- opt$options$gwas_code
chr <- opt$options$chr
ref_path <- opt$options$ref_path
n <- opt$options$n_snps

print(paste0("reading gwas: ", key, " - chr: ", chr))

ref_panel = fread(paste0(ref_path,"/",chr,".bim"))
names(ref_panel) =c( "chr","rsID", "MAF", "pos", "A1", 'A2')

typed = fread(paste0("z_scores/z_",key,"_",chr,".txt"))

setkey(ref_panel, 'rsID')
setkey(typed, 'rsID')

# add MAF to typed SNPs
typed$MAF = ref_panel[typed, MAF]

# calculate actual MAF
get_maf <- function(x){min(x['MAF'], 1-x['MAF'])}
typed$MAF_minor <- apply(typed, 1, get_maf)

Q_MAF = c(0, 0.01, 0.05, 0.1, 0.5, 1)

typed$MAF_bin <- cut(typed$MAF_minor, Q_MAF, labels = F)

n_per_bin <- n / (length(Q_MAF)-1 )

sampled <- typed %>% 
  group_by(MAF_bin) %>% 
  sample_n(n_per_bin) %>%
  select(rsID, pos, A0, A1, Z, P)

fwrite(sampled, paste0("z_scores/z_",key,"_",chr,"_sampled.txt"), quote = F, sep = '\t')