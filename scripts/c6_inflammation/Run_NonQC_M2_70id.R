##Model 2
#Run_NonQC_M2_70id.R

rm(list=ls())
library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)

meth <- fread("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/maggie/NewID/beta_enmix_all_twins_blood_450k.csv", data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

info <- read.csv("/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/maggie/NewID/15012020_New_5y_70id_info_fib.csv")
info[12] <- "F"
meth <- meth[, match(info$Barcode, colnames(meth))]

info$Zygosity[info$Zygosity == "MZ"] <- info$KCLfam[info$Zygosity == "MZ"]
info$Zygosity[info$Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Zygosity == "DZ"])

info$Fibrinogen_log <- as.numeric(info$Fibrinogen_log)
info$Age <- as.numeric(info$Age)
info$BMI_Nearest <- as.numeric(info$BMI_Nearest)
info$Smoking_Nearest2 <- as.factor(info$Smoking_Nearest2)
info$CD8T <- as.numeric(info$CD8T)
info$CD4T <- as.numeric(info$CD4T)
info$NK <- as.numeric(info$NK)
info$Bcell <- as.numeric(info$Bcell)
info$Mono <- as.numeric(info$Mono)
info$Gran <- as.numeric(info$Gran)
info$Sentrix_ID <- as.factor(info$Sentrix_ID)
info$Sentrix_Position <- as.factor(info$Sentrix_Position)
info$KCLfam <- as.factor(info$KCLfam)
info$Zygosity <- as.factor(info$Zygosity)
info$date_diff <- as.numeric(info$date_diff)

f1 <- y ~ Fibrinogen_log + Age + BMI_Nearest + Smoking_Nearest2 + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|KCLfam) + (1|Zygosity) + date_diff         
f2 <- y ~ Age + BMI_Nearest + Smoking_Nearest2 + CD8T + CD4T + NK + Bcell + Mono + Gran + (1|Sentrix_ID) + (1|Sentrix_Position) + (1|KCLfam) + (1|Zygosity) + date_diff

ewasOne <- function(x){
  info$y = as.numeric(meth[x,]) 
  fit = lmer(f1, data = info) 
  null = lmer(f2, data = info) 
  a = anova(fit,null)
  tbl = summary(fit)$coef 
  result = c(tbl[2, 1], tbl[2, 2], tbl[2, 3], a$"Pr(>Chisq)"[2], summary(fit)$devcomp$dims[1])
}

library(parallel) 
options(mc.cores = 4) 
resu <- do.call(rbind, mclapply(1:NROW(meth), ewasOne))
resu_out <- cbind(rownames(meth),resu) 
colnames(resu_out) <- c("CpG","Beta","SE","T-value","P-value","N")
resu_out <- as.data.frame(resu_out)

write.csv(resu_out, "/mnt/lustre/groups/bell/epigenetics/Analysis/subprojects/maggie/NewID/NonQC_M2_70id_fib_Result.csv", quote = F, row.names = F)

