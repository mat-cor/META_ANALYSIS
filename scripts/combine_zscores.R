#!/usr/bin/env Rscript

# packs <- c("optparse","data.table","R.utils","dplyr")
# for (p in packs) {
#   if( !require(p, character.only = T)) {
#     print(p)
#     install.packages( p,  repos = c(CRAN = "https://cloud.r-project.org") )
#     require(p, character.only = T)
#   }
# }

require(optparse)
require(data.table)
require(R.utils)
require(dplyr)

option_list = list(
  make_option(c("-s", "--original_sstat"), type="character", default=NULL,
              help="original sumstat", metavar="character"),
  make_option(c("-i", "--imputed_sstat"), type="character", default=NULL,
              help="imputed sumstat", metavar="character"),
  make_option(c("-o", "--out"), type="character", default = NULL,
              help="output file name ", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

o <- opt$options$original_sstat
print(paste("reading file:", o))
ori <- fread(o, header=T)

i <- opt$options$imputed_sstat
print(paste("reading file:", i))
imp <- fread(i, header=T)

colnames(ori)[1] <- 'CHR'
ori <- ori %>%
  mutate(Z = BETA/SE,
         rsID = paste0(CHR, ':', POS, ':', Allele1, ':', Allele2),
         raiss.imputed = 0) %>%
  select(rsID, '#CHR' = CHR, POS, Allele1, Allele2, Z, p.value, raiss.imputed, BETA, SE, AF_Allele2, imputationInfo, N)

head(ori)

imp <- imp %>%
  filter(!rsID %in% ori$rsID) %>%
  mutate('#CHR' = as.integer(sub(":.*", "", rsID)),
         p.value = 2*pnorm(-abs(Z)),
         raiss.imputed = 1) %>%
  select(rsID, '#CHR', POS = pos, Allele1 = A0, Allele2 = A1, Z, p.value, raiss.imputed, Var, ld_score)

head(imp)

out <- bind_rows(ori, imp)
out <- out %>%
  arrange(`#CHR`, POS)

oo <- opt$options$out
fwrite(out, oo, sep = '\t')