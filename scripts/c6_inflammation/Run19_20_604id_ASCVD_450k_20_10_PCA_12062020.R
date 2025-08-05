rm(list=ls())
library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)

meth <- fread("/scratch/groups/bell/epigenetics/Analysis/subprojects/maggie/450k/450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv", data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

info <- read.csv("/scratch/groups/bell/epigenetics/Analysis/subprojects/604id_info_20PCA_12062020.csv")
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
info$PC1 <- as.numeric(info$PC1)
info$PC2 <- as.numeric(info$PC2)
info$PC3 <- as.numeric(info$PC3)
info$PC4 <- as.numeric(info$PC4)
info$PC5 <- as.numeric(info$PC5)
info$PC6 <- as.numeric(info$PC6)
info$PC7 <- as.numeric(info$PC7)
info$PC8 <- as.numeric(info$PC8)
info$PC9 <- as.numeric(info$PC9)
info$PC10 <- as.numeric(info$PC10)
info$PC11 <- as.numeric(info$PC11)
info$PC12 <- as.numeric(info$PC12)
info$PC13 <- as.numeric(info$PC13)
info$PC14 <- as.numeric(info$PC14)
info$PC15 <- as.numeric(info$PC15)
info$PC16 <- as.numeric(info$PC16)
info$PC17 <- as.numeric(info$PC17)
info$PC18 <- as.numeric(info$PC18)
info$PC19 <- as.numeric(info$PC19)
info$PC20 <- as.numeric(info$PC20)

f1 <- y ~ ASCVD_100 + Age + BMI_Nearest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + (1|KCLfam)     
f2 <- y ~ Age + BMI_Nearest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + (1|KCLfam)  

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
colnames(resu_out) <- c("CpG","Beta","SE","T-value","P-value","N")
resu_out <- as.data.frame(resu_out)

write.csv(resu_out, "/scratch/groups/bell/epigenetics/Analysis/subprojects/maggie/ASCVD/results/19_604id_ASCVD_450k_20PCA_12062020.csv", quote = F, row.names = F)

#####
rm(list=ls())
library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)

meth <- fread("/scratch/groups/bell/epigenetics/Analysis/subprojects/maggie/450k/450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv", data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

info <- read.csv("/scratch/groups/bell/epigenetics/Analysis/subprojects/604id_info_20PCA_12062020.csv")
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
info$PC1 <- as.numeric(info$PC1)
info$PC2 <- as.numeric(info$PC2)
info$PC3 <- as.numeric(info$PC3)
info$PC4 <- as.numeric(info$PC4)
info$PC5 <- as.numeric(info$PC5)
info$PC6 <- as.numeric(info$PC6)
info$PC7 <- as.numeric(info$PC7)
info$PC8 <- as.numeric(info$PC8)
info$PC9 <- as.numeric(info$PC9)
info$PC10 <- as.numeric(info$PC10)
info$PC11 <- as.numeric(info$PC11)
info$PC12 <- as.numeric(info$PC12)
info$PC13 <- as.numeric(info$PC13)
info$PC14 <- as.numeric(info$PC14)
info$PC15 <- as.numeric(info$PC15)
info$PC16 <- as.numeric(info$PC16)
info$PC17 <- as.numeric(info$PC17)
info$PC18 <- as.numeric(info$PC18)
info$PC19 <- as.numeric(info$PC19)
info$PC20 <- as.numeric(info$PC20)

f1 <- y ~ ASCVD_100 + Age + BMI_Nearest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1|KCLfam)     
f2 <- y ~ Age + BMI_Nearest + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (1|KCLfam)  

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
colnames(resu_out) <- c("CpG","Beta","SE","T-value","P-value","N")
resu_out <- as.data.frame(resu_out)

write.csv(resu_out, "/scratch/groups/bell/epigenetics/Analysis/subprojects/maggie/ASCVD/results/20_604id_ASCVD_450k_10PCA_12062020.csv", quote = F, row.names = F)