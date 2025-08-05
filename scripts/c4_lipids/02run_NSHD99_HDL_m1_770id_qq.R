rm(list=ls())
library(lme4)
library(data.table)
library(dplyr)
options(stringsAsFactors=FALSE)

meth <- fread("/scratch/users/k1787261/NSHD/NSHD_99_beta_enmix_BlExcluded_SexChrExcluded_750119cg_1460basename.csv", data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

info <- read.csv("/scratch/users/k1787261/rory_lipids/info/99_NSHD_info_lipid_770id.csv")
meth <- meth[, match(info$Basename, colnames(meth))]

info$HDLmgdl <- as.numeric(info$HDLmgdl)
info$BMI <- as.numeric(info$BMI)
info$smok <- as.factor(info$smok)
info$CD8T <- as.numeric(info$CD8T)
info$CD4T <- as.numeric(info$CD4T)
info$NK <- as.numeric(info$NK)
info$Bcell <- as.numeric(info$Bcell)
info$Mono <- as.numeric(info$Mono)
info$Gran <- as.numeric(info$Gran)
info$Sentrix_ID <- as.factor(info$Sentrix_ID)
info$Sentrix_Position <- as.factor(info$Sentrix_Position)
info$sex <- as.factor(info$sex)
info$medchol_53x <- as.factor(info$medchol_53x)

f1 <- y ~ HDLmgdl + sex + BMI + smok + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + medchol_53x
f2 <- y ~ sex + BMI + smok + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + medchol_53x

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
options(mc.cores = 48) 
resu <- do.call(rbind, mclapply(1:NROW(meth), ewasOne))
resu_out <- cbind(rownames(meth), resu) 
colnames(resu_out) <- c("cpg","beta","se","t_value","p_val","N")
resu_out <- as.data.frame(resu_out)

write.csv(resu_out, "/scratch/users/k1787261/rory_lipids/results/qq/qq102_NSHD99_HDL_medco_770id.csv", quote = F, row.names = F)
