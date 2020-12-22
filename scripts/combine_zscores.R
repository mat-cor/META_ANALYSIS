#!/usr/bin/env Rscript

packs <- c("optparse","data.table","R.utils","dplyr")

for (p in packs) {
  if( !require(p, character.only = T)) {
    print(p)
    install.packages( p,  repos = c(CRAN = "https://cloud.r-project.org") )  
  }
}


option_list = list(
  make_option(c("-ori", "--original-sstat"), type="character", default=NULL,
              help="original sumstat", metavar="character"),
  make_option(c("-imp", "--imputed-sstat"), type="character", default=NULL,
              help="original sumstat", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

o <- 'Projects/covid19-hgi/test_imputation/FinnGen.Karjalainen.ANA_C2_V2.11.ALL.ALL.FIN.357.238354.SAIGE.20200915.txt.munged.AF.0.001.INFO.0.6.gz.gnomad_v3_b38_ref_fin.gz.gnomad_v2.1.1_b37_ref_nfe.gz'
o <- opt$options$ori
print(paste("reading file:", o))
ori <- fread(o, header=T)

i <- 'Projects/covid19-hgi/test_imputation/FinnGen.Karjalainen.ANA_C2_V2.11.ALL.ALL.FIN.357.238354.SAIGE.20200915.txt.munged.AF.0.001.INFO.0.6.gz.gnomad_v3_b38_ref_fin.gz.imputed.txt.gz'
i <- opt$options$imp
print(paste("reading file:", i))
imp <- fread(i, header=T)

colnames(ori)[1] <- 'CHR'
ori <- ori %>%
  mutate(Z = BETA/SE,
         rsID = paste0(CHR, ':', POS, ':', Allele1, ':', Allele2),
         raiss.imputed = 0) %>%
  select(rsID, '#CHR' = CHR, POS, Allele1, Allele2, Z, p.value, raiss.imputed, BETA, SE, AF_gnomad_v2.1.1_b37_ref_nfe, AF_fc, imputationInfo, N)

head(ori)

imp <- imp %>%
  filter(!rsID %in% ori$rsID) %>%
  mutate('#CHR' = as.integer(sub(":.*", "", rsID)),
         p.value = 2*pnorm(-abs(Z)),
         raiss.imputed = 1) %>%
  select(rsID, '#CHR', POS = pos, Allele1 = A0, Allele2 = A1, Z, p.value, Var, ld_score)

head(imp)

out <- bind_rows(ori, imp)

is.data.frame(out)

out <- out %>%
  arrange(`#CHR`, POS)

oo <- opt$options$o
fwrite(out, oo, sep = '\t')