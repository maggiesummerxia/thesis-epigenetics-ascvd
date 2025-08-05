rm(list=ls())
library(lme4)
library(tidyverse)
options(stringsAsFactors=FALSE)

meth <- read_rds("/scratch/users/k1787261/EPIC/beta_norm_qced_ncds_new_1359_blrm_sexrm.rds")
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

# /Users/maggiexia/OneDrive/maggie_phd/ncds_qc/data/processed/beta_norm_qced_ncds_new_1359_blrm_sexrm.rds
# /scratch/users/k1787261/EPIC/beta_norm_qced_ncds_new_1359_blrm_sexrm.rds

info <- read.csv("/scratch/users/k1787261/meno/info/ncds_new_meno_status_info_417id.csv")
meth <- meth[, match(info$basename, colnames(meth))]

# /Users/maggiexia/OneDrive/maggie_phd/ncds_qc/data/processed/ncds_new_meno_status_info_417id.csv
# /scratch/users/k1787261/meno/info/ncds_new_meno_status_info_417id.csv

info$chip_id <- as.factor(info$chip_id)
info$smoking <- as.factor(info$smoking50)
info$meno_status <- as.factor(info$meno_status)

f1 <- y ~ meno_status + smoking + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|chip_id) + (1|chip_position)
f2 <- y ~ smoking + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|chip_id) + (1|chip_position)

ewasOne <- function(x){
  info$y = as.numeric(meth[x,]) 
  info$y = qqnorm(info$y, plot.it=FALSE)$x
  fit = lmer(f1, data = info) 
  null = lmer(f2, data = info) 
  a = anova(fit,null)
  tbl = summary(fit)$coef 
  result = c(tbl[2, 1], tbl[2, 2], tbl[2, 3], a$"Pr(>Chisq)"[2], summary(fit)$devcomp$dims[1])
}

library(parallel) 
options(mc.cores = 12) 
resu <- do.call(rbind, mclapply(1:NROW(meth), ewasOne))
resu_out <- cbind(rownames(meth), resu) 
colnames(resu_out) <- c("MARKER","BETA","SE","TSTATISTIC","PVAL","N")
resu_out <- as.data.frame(resu_out)
resu_out$PVAL <- as.numeric(resu_out$PVAL)
resu_out$FDR <- p.adjust(resu_out$PVAL, method = "fdr")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
Locations <- as.data.frame(Locations)
Locations$MARKER <- rownames(Locations)
Locations <- Locations[, c("MARKER", "chr", "pos")]
colnames(Locations) <- c("MARKER", "CHR", "POS")

resu_out <- left_join(resu_out, Locations, by = "MARKER")
resu_out <- resu_out[, c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")]

write_tsv(resu_out, "/scratch/users/k1787261/meno/results/ncds_new_meno_status_resu_417id.txt")


