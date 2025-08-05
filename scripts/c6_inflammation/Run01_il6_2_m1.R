rm(list=ls())
library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)

meth <- fread("/scratch/users/k1787261/450k/450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv", data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

info <- read.csv("/scratch/users/k1787261/cytokines/info/info/239id_il6_2_info_450k_blood.csv")
meth <- meth[, match(info$Barcode, colnames(meth))]

info$Zygosity[info$Zygosity == "MZ"] <- info$KCLfam[info$Zygosity == "MZ"]
info$Zygosity[info$Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Zygosity == "DZ"])

info$il6_2 <- as.numeric(info$il6_2)
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
info$Date_diff <- as.numeric(info$Date_diff)

f1 <- y ~ il6_2 + Age + BMI_Nearest + Smoking_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|Zygosity) + (1|KCLfam) + Date_diff
f2 <- y ~ Age + BMI_Nearest + Smoking_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|Zygosity) + (1|KCLfam) + Date_diff

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
colnames(resu_out) <- c("CpG","BETA","SE","t_value","p_value","N")
resu_out <- as.data.frame(resu_out)
resu_out$FDR <- p.adjust(resu_out$P.value, method = "fdr")
resu_out <- arragne(resu_out, p.value)

write.csv(resu_out, "/scratch/users/k1787261/cytokines/results/resu01_il6_2_239id_450k_blood_m1.csv", quote = F, row.names = F)

