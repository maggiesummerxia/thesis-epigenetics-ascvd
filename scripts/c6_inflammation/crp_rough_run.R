library(data.table)
library(here)
library(tidyverse)
library(readxl)
library(lme4)

meth <- fread(here("master_data", "datasets", "TwinsUK", "450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv"), data.table = FALSE)
info <- read.csv(here("crp", "data", "info_b1andb2_crp.csv"))
# crp_sig <- read_excel(here("crp", "data", "41467_2022_29792_MOESM4_ESM.xlsx"), sheet = 3)
# crp_sig <- crp_sig[, 1]

# meth <- meth %>%
#  inner_join(crp_sig, by = c("cg" = "ID"))

rownames(meth) <- meth[, 1]
meth <- meth[, -1]

meth_id <- tibble("Basename" = colnames(meth))
crp_id <- tibble("Basename" = info$Basename)
ol_id <- inner_join(meth_id, crp_id, by = "Basename") #415
info <- info %>%
inner_join(ol_id, by = "Basename")

meth <- meth[, match(info$Basename, colnames(meth))]

info$Actual_Zygosity[info$Actual_Zygosity == "MZ"] <- info$KCLfam[info$Actual_Zygosity == "MZ"]
info$Actual_Zygosity[info$Actual_Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Actual_Zygosity == "DZ"])

info$newsmoke <- as.factor(info$newsmoke)
info$plate <- as.factor(info$plate)
info$position <- as.factor(info$position)
info$Actual_Zygosity <- as.factor(info$Actual_Zygosity)
info$KCLfam <- as.factor(info$KCLfam)
info$Batch1 <- as.factor(info$Batch1)

f1 <- y ~ logcrp + newsmoke + Batch1 + (1|KCLfam) + (1|plate) + (1|position) + 
  BMI + Age + (1|ACTUAL_ZYGOSITY) + CD8T + CD4T + NK + Bcell + Mono + Gran

f2 <- y ~ newsmoke + Batch1 + (1|KCLfam) + (1|plate) + (1|position) + 
  BMI + Age + (1|ACTUAL_ZYGOSITY) + CD8T + CD4T + NK + Bcell + Mono + Gran

ewasOne <- function(x){
  info$y = as.numeric(meth[x,]) 
  fit = lmer(f1, data = info) 
  null = lmer(f2, data = info) 
  a = anova(fit,null)
  tbl = summary(fit)$coef 
  result = c(tbl[2, 1], tbl[2, 2], tbl[2, 3], a$"Pr(>Chisq)"[2], summary(fit)$devcomp$dims[1])
}

library(parallel) 
options(mc.cores = 6) 
resu <- do.call(rbind, mclapply(1:NROW(meth), ewasOne))
resu_out <- cbind(rownames(meth), resu) 
colnames(resu_out) <- c("CpG","beta","SE", "t_value", "p","N")
resu_out <- as.data.frame(resu_out)
resu_out$p <- as.numeric(resu_out$p)
resu_out <- resu_out %>%
  arrange(p)

write_csv(resu_out, here("crp", "results", "resu_crp_tuk_415id.csv"))

# explore date diff
info$date_diff <- as.Date(info$DATE_OF_VISIT, format = "%d/%m/%Y") - as.Date(info$DNAextraction, format = "%d/%m/%Y")
# 10 year