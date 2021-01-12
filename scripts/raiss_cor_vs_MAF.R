require(optparse)
require(data.table)
require(R.utils)
require(ggplot2)

option_list = list(
  make_option(c("-g", "--gwas_code"), type="character", default=NULL,
              help="gwas code", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL,
              help="chromosome", metavar="character"),
  make_option(c("-r", "--ref_path"), type="character", default=NULL,
              help="ref panel path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

key <- opt$options$gwas_code
chr <- opt$options$chr
ref_path <- opt$options$ref_path

print(paste0("reading gwas: ", key, " - chr: ", chr))

ref_panel = fread(paste0(ref_path,chr,".bim"))
typed = fread(paste0("./z_scores/z_",key,"_",chr,".txt"))

masked = fread(paste0("./z_scores_masked/z_scores_masked_z_",key,"_",chr,".txt"))

imputed = fread(paste0("./z_scores_imputed_test/z_",key,"_",chr,"_0.0001.txt"))

dim(typed)
dim(masked)
dim(imputed)

names(ref_panel) =c( "chr","rsID", "MAF", "pos", "A1", 'A2')

masked_SNP = setdiff(typed$rsID, masked$rsID)
length(masked_SNP)
setkey(ref_panel, 'rsID')
setkey(typed, 'rsID')
setkey(masked, 'rsID')
setkey(imputed, 'rsID')

head(ref_panel[masked_SNP, MAF ])
ggplot(ref_panel, aes(x=MAF)) + geom_histogram()

perf_table = data.frame(Z_typed = typed[masked_SNP, Z], Z_imputed=imputed[masked_SNP, Z], R2=1-imputed[masked_SNP, Var], MAF = ref_panel[masked_SNP, MAF ], nchar = ref_panel[masked_SNP, nchar(paste0(A1, A2)) ])
perf_table = perf_table[!is.na(perf_table$Z_imputed),]
head(perf_table)
cor(perf_table[perf_table$nchar==2,1:2])
cor(perf_table[,1:2])

get_maf <- function(x){min(x['MAF'], 1-x['MAF'])}

perf_table["MAF_minor"] = apply(perf_table, 1, get_maf)

Q_MAF = c(0, 0.01, 0.05, 0.1, 0.5) # quantile(perf_table$MAF_minor, p= seq(0,1, 0.25))


perf_table["MAF_quantile"] = ""
for(j in 2:length(Q_MAF)){
  print(cor(perf_table[ (perf_table$MAF_minor > Q_MAF[j-1]) & (perf_table$MAF_minor < Q_MAF[j]),1:2]))
  cpi = round((cor(perf_table[ (perf_table$MAF_minor > Q_MAF[j-1]) & (perf_table$MAF_minor < Q_MAF[j]),1:2])[1,2]),2)
  print(cpi)
  perf_table[((perf_table$MAF_minor >= Q_MAF[j-1]) & (perf_table$MAF_minor <= Q_MAF[j])), "MAF_quantile"] = paste(round(Q_MAF[j-1], 2)*100,"% < MAF <", round(Q_MAF[j], 2)*100, "%")
}
head(perf_table)
R2_bins = c(0.6,0.8, 1.0)
perf_table["R2_quantile"] = ""
for(j in 2:length(R2_bins)){
  perf_table[((perf_table$R2 >= R2_bins[j-1]) & (perf_table$R2 <= R2_bins[j])), "R2_quantile"] = paste(round(R2_bins[j-1], 2),"< R2 <", round(R2_bins[j], 2))
}
head(perf_table)


perf_table["L1_error"] = abs(perf_table$Z_typed-perf_table$Z_imputed)
#unique(perf_table$MAF_quantile)
#ordered(MAF_quantile, levels=c('0 % < MAF < 1 %', '1 % < MAF < 5 %','5 % < MAF < 10 %', '10 % < MAF < 100 %'))
p = ggplot(perf_table, aes(x=  ordered(MAF_quantile, levels=c('0 % < MAF < 1 %', '1 % < MAF < 5 %','5 % < MAF < 10 %', '10 % < MAF < 50 %')), fill=R2_quantile, y=L1_error)) + geom_boxplot( outlier.shape = NA, notch=TRUE)
p = p + ylim(c(0,0.75)) + coord_flip() + scale_fill_manual(values=c( "orange", "royalblue"))
p = p + theme(legend.position="top") + xlab("")
ggsave(paste0("./diag/L1_error_",key,".png"), width=6, height=6,plot=p)

p = ggplot(perf_table, aes(x=  ordered(MAF_quantile, levels=c('0 % < MAF < 1 %', '1 % < MAF < 5 %','5 % < MAF < 10 %', '10 % < MAF < 50 %')), fill=R2_quantile, y=L1_error))+  geom_violin()
p = p + ylim(c(0,0.75)) + coord_flip()
p = p + theme(legend.position="top") + xlab("")
p + stat_summary(fun.data="mean_sdl", mult=1,
                 geom="crossbar", width=0.2 )

ggsave(paste0("./diag/L1_error_",key,"violin.png"), width=6, height=6,plot=p)


present_maf = intersect(c('0 % < MAF < 1 %', '1 % < MAF < 5 %','5 % < MAF < 10 %', '10 % < MAF < 50 %'),unique(perf_table$MAF_quantile))
nbins = length(present_maf)

print(present_maf)
Cor_tab = data.frame(row.names=present_maf, cor=rep(0, nbins), conf_interval= rep(0, nbins))

for(maf_b in present_maf){
  ct = cor.test(perf_table[which(perf_table$MAF_quantile == maf_b), "Z_typed" ], perf_table[which(perf_table$MAF_quantile == maf_b), "Z_imputed"])
  
  Cor_tab[maf_b, "cor"] = signif(ct$estimate, 2)
  Cor_tab[maf_b, "conf_interval"] = signif(ct$conf.int[2] - ct$estimate, 2)
}
write.table(Cor_tab, paste0("./diag/cor_",key,".txt"))
