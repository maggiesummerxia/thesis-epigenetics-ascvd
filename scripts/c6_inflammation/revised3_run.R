meth_name <- "/scratch/prj/bell_epigen/maggie/master_data/450k_blood_all_twins_beta_enmix_probeqced_990ids_ready_sex_chr_excluded.csv"

info_name <- "/scratch/prj/bell_epigen/maggie/cytoadipokine/data/il6_rerun_info_191id_LS.csv"

f1 <- y ~ il6_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff

f2 <- y ~ Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff

resu_name <- "/scratch/prj/bell_epigen/maggie/cytoadipokine/results/revised3_il6_base_twinsuk_20221012.csv"

source("/scratch/prj/bell_epigen/maggie/cytoadipokine/scripts/master_script.R")

# il6_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + Smoking_Nearest

# il6_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + crp_mgl

# leptin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff

# leptin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + Smoking_Nearest

# leptin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + BMI_Nearest

# adiponectin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff

# adiponectin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + Smoking_Nearest

# adiponectin_log + Age + CD8T + CD4T + NK + Bcell + Mono + (1|Sentrix_ID) + (1|Zygosity) + (1|KCLfam) + Date_diff + BMI_Nearest


library(lme4)
library(data.table)
options(stringsAsFactors=FALSE)

meth <- fread(meth_name, data.table=FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]

#test
#meth <- meth[1:100, ]
#

info <- read.csv(info_name)
meth <- meth[, match(info$Barcode, colnames(meth))]

info$Zygosity[info$Zygosity == "MZ"] <- info$KCLfam[info$Zygosity == "MZ"]
info$Zygosity[info$Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Zygosity == "DZ"])

info$Smoking_Nearest <- as.factor(info$Smoking_Nearest)
info$Sentrix_ID <- as.factor(info$Sentrix_ID)
info$Sentrix_Position <- as.factor(info$Sentrix_Position)
info$Zygosity <- as.factor(info$Zygosity)
info$KCLfam <- as.factor(info$KCLfam)

ewasOne <- function(x){
  info$y = as.numeric(meth[x,]) 
  fit = lmer(f1, data = info) 
  null = lmer(f2, data = info) 
  a = anova(fit,null)
  tbl = summary(fit)$coef 
  mean_meth <- mean(info$y, na.rm = TRUE)
  sd_meth <- sd(info$y, na.rm = TRUE)
  result = c(tbl[2, 1], tbl[2, 2], tbl[2, 3], a$"Pr(>Chisq)"[2], summary(fit)$devcomp$dims[1], mean_meth, sd_meth)
}

library(parallel) 
options(mc.cores = 48) 
resu <- do.call(rbind, mclapply(1:NROW(meth), ewasOne))
resu_out <- cbind(rownames(meth), resu) 
colnames(resu_out) <- c("CpG","beta","SE", "t_value", "p","N", "meth_beta_mean", "meth_beta_sd")
resu_out <- as.data.frame(resu_out)

resu_out <- resu_out[, c("CpG", "beta", "SE", "p", "N", "meth_beta_mean", "meth_beta_sd")]

write.csv(resu_out, resu_name, quote = F, row.names = F)

rm(list=ls())