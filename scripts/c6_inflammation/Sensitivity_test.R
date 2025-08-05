rm(list=ls())
library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)
library(readr)
library(dplyr)
info <- read_csv('717id_info_450k_ASCVD_5y.csv')
#info <- info %>% filter(ASCVD_100 < 30)
meth <- fread('../../../Documents/Rosa_trans/450k_beta_subset_topcpg_21052020.csv', data.table = FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]
meth <- meth[, match(info$Barcode, colnames(meth))]

info$Zygosity[info$Zygosity == "MZ"] <- info$KCLfam[info$Zygosity == "MZ"]
info$Zygosity[info$Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Zygosity == "DZ"])

info$ASCVD_100 <- as.numeric(info$ASCVD_100)
info$Age <- as.numeric(info$Age)
info$BMI_Nearest <- as.numeric(info$BMI_Nearest)
info$Smoking_Nearest <- as.factor(info$Smoking_Nearest)
info$CD8T <- as.numeric(info$CD8T)
info$CD4T <- as.numeric(info$CD4T)
info$NK <- as.numeric(info$NK)
info$Bcell <- as.numeric(info$Bcell)
info$Mono <- as.numeric(info$Mono)
info$Gran <- as.numeric(info$Gran)
info$Sentrix_ID <- as.factor(info$Sentrix_ID)
info$Sentrix_Position <- as.factor(info$Sentrix_Position)
info$Zygosity <- as.factor(info$Zygosity)
info$KCLfam <- as.factor(info$KCLfam)

f1 <- y ~ ASCVD_100 + Age + Smoking_Nearest + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|Zygosity) + (1|KCLfam)
f2 <- y ~ Age + Smoking_Nearest + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|Zygosity) + (1|KCLfam)

ewasOne <- function(x){
  info$y = as.numeric(meth[x,]) 
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
colnames(resu_out) <- c("CpG","Beta","SE","T-value","P-value","N")
resu_out <- as.data.frame(resu_out)
#write_csv(resu_out, 'Runtop5_score_less_than_30.csv')

resu_out_5 <- resu_out %>% filter(CpG == "cg26699283" | CpG == "cg01773518" | CpG == "cg04925511" | CpG == "cg18696321" | CpG == "cg20241658")
resu_out_5$Beta <- round(as.numeric(resu_out_5$Beta), 3)
resu_out_5 <- resu_out_5 %>% arrange(`P-value`)
