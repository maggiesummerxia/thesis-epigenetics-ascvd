library(here)
library(tidyverse)
library(data.table)
#library(pwr)

# power simulation function ----
power_simulation <- function(case_n, control_n, rep_time, 
                             case_df = info_case, control_df = info_control,
                             full_df = info,
                             barcode_colname = "Basename", cg_number,
                             model_norm, model_raw){
  
  result_simu <- matrix(, nrow = rep_time, ncol = 10)
  
  x = 1
  
  repeat {
    case <- sample(case_df[[barcode_colname]], case_n, replace = TRUE, prob = NULL)
    control <- sample(control_df[[barcode_colname]], control_n, replace = TRUE, prob = NULL)
    
    all <- c(case, control)
    info_simulation <- data.frame("Barcode" = all)
    info_simulation <- left_join(info_simulation, full_df, by = c("Barcode" = barcode_colname))
    info_simulation$norm <- qqnorm(info_simulation[[cg_number]], plot.it=FALSE)$x
    
    full <- lm(model_norm, data = info_simulation)
    n <- nobs(full)
    fitted <- summary(full)$coefficients
    
    full_raw <- lm(model_raw, data = info_simulation)
    fitted_raw <- summary(full_raw)$coefficients
    
    mean_case_raw <- mean(info_simulation[info_simulation$ASCVD == "1", ][[cg_number]], na.rm = TRUE)
    mean_control_raw <- mean(info_simulation[info_simulation$ASCVD == "0", ][[cg_number]], na.rm = TRUE)
    mean_diff_raw <- mean_case_raw - mean_control_raw
    
    results <- c(fitted[2,1],fitted[2,2],fitted[2,4],fitted_raw[2,1],fitted_raw[2,2],fitted_raw[2,4],mean_case_raw,mean_control_raw,mean_diff_raw,n)
    result_simu[x, ] <- results
    colnames(result_simu) <- c("beta","se","p","beta_raw","se_raw","p_raw","mean_case_raw","mean_control_raw","mean_diff_raw","n")
    result_simu <- data.frame(result_simu)
    
    if (x > rep_time -1){
      break
    }
    
    x = x + 1
  }
  
  result_simu
}

# twinsuk prep ----
info <- read_csv(here("ichd", "info", "TwinsUK_450k_iCVD_67cases_535controls.csv"))

info$Zygosity[info$Zygosity == "MZ"] <- info$KCLfam[info$Zygosity == "MZ"]
info$Zygosity[info$Zygosity == "DZ"] <- paste0("ID", info$KCLid[info$Zygosity == "DZ"])

factor_cols <- c("Sex", "Smoking_Nearest", "Sentrix_ID", "Sentrix_Position", "Zygosity", "KCLfam", "ASCVD", "Status")
info[factor_cols] <- lapply(info[factor_cols], factor)

meth <- fread(here("master_data", "datasets", "TwinsUK", "450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv"), data.table = FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]
meth <- meth[, match(info$Barcode, colnames(meth))]
meth_cg04580929 <- meth %>% filter(rownames(meth) == "cg04580929")
meth_cg04580929 <- data.frame(t(meth_cg04580929))
meth_cg09117673 <- meth %>% filter(rownames(meth) == "cg09117673")
meth_cg09117673 <- data.frame(t(meth_cg09117673))

info$cg04580929 <- meth_cg04580929$cg04580929
info$cg09117673 <- meth_cg09117673$cg09117673


# cg04580929	-1.678	0.338	7.04E-07	1	601	-??	5	PDE4D	450k only
# cg09117673	0.702	0.144	1.05E-06	3	1683	+++	7	WBSCR28	450k & EPIC

ggplot(info, aes(ASCVD, cg04580929)) +
  geom_boxplot()

ggplot(info, aes(ASCVD, cg09117673)) +
  geom_boxplot()

# check qq normalised values
info$cg04580929_norm <- qqnorm(info$cg04580929, plot.it=FALSE)$x
info$cg09117673_norm <- qqnorm(info$cg09117673, plot.it=FALSE)$x

ggplot(info, aes(ASCVD, cg04580929_norm)) +
  geom_boxplot()
ggplot(info, aes(ASCVD, cg09117673_norm)) +
  geom_boxplot()

# check if it returns the same result as previous analyses
full <- lm(cg04580929_norm ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam, data = info)
n <- nobs(full)
fitted <- summary(full)$coefficients
results=c(fitted[2,1],fitted[2,2],fitted[2,4],n)

# power analysis dataframes
info_case <- info %>% filter(ASCVD == "1")
info_control <- info %>% filter(ASCVD == "0")

# twinsuk cg04580929 ----
set.seed(1701)

# cg04580929
cg04580929_simu_67_535_1000 <- power_simulation(case_n = 67, control_n = 535, rep_time = 1000,
                                     case_df = info_case, control_df = info_control,
                                     full_df = info,barcode_colname = "Barcode", cg_number = "cg04580929",
                                     model_norm = norm ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam, 
                                     model_raw = cg04580929 ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam)

cg04580929_simu_67_535_1000 <- cg04580929_simu_67_535_1000 %>%
  mutate(
    sig = case_when(
      p < 1.186969e-07 ~ 1,
      p >= 1.186969e-07 ~ 0
    )
  )

table(cg04580929_simu_67_535_1000$sig)

write_csv(cg04580929_simu_67_535_1000, here("ichd_followup", "results", "power", "cg04580929_tuk_67_535_1000.csv"))

# twinsuk cg09117673 ----

# cg09117673
set.seed(1701)

cg09117673_simu_67_535_1000 <- power_simulation(case_n = 67, control_n = 535, rep_time = 1000,
                                                case_df = info_case, control_df = info_control,
                                                full_df = info,barcode_colname = "Barcode", cg_number = "cg09117673",
                                                model_norm = norm ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam, 
                                                model_raw = cg09117673 ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam)

cg09117673_simu_67_535_1000 <- cg09117673_simu_67_535_1000 %>%
  mutate(
    sig = case_when(
      p < 1.186969e-07 ~ 1,
      p >= 1.186969e-07 ~ 0
    )
  )

table(cg09117673_simu_67_535_1000$sig)

#  0   1 
# 780 220 


set.seed(1701)

cg09117673_simu_134_1070_1000 <- power_simulation(case_n = 134, control_n = 1070, rep_time = 1000,
                                                case_df = info_case, control_df = info_control,
                                                full_df = info,barcode_colname = "Barcode", cg_number = "cg09117673",
                                                model_norm = norm ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam, 
                                                model_raw = cg09117673 ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam)

cg09117673_simu_134_1070_1000 <- cg09117673_simu_134_1070_1000 %>%
  mutate(
    sig = case_when(
      p < 1.186969e-07 ~ 1,
      p >= 1.186969e-07 ~ 0
    )
  )

table(cg09117673_simu_134_1070_1000$sig)

# 0   1 
# 551 449 

write_csv(cg09117673_simu_134_1070_1000, here("ichd_followup", "results", "power", "cg09117673_tuk_134_1070_1000.csv"))

set.seed(1701)

cg09117673_simu_268_2140_1000 <- power_simulation(case_n = 268, control_n = 2140, rep_time = 1000,
                                                  case_df = info_case, control_df = info_control,
                                                  full_df = info,barcode_colname = "Barcode", cg_number = "cg09117673",
                                                  model_norm = norm ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam, 
                                                  model_raw = cg09117673 ~ ASCVD + Age + BMI_Nearest + CD8T + CD4T + NK + Bcell + Mono + Gran + Sex + Smoking_Nearest + Status + FollowUp + Sentrix_ID + Sentrix_Position + Zygosity + KCLfam)

cg09117673_simu_268_2140_1000 <- cg09117673_simu_268_2140_1000 %>%
  mutate(
    sig = case_when(
      p < 1.186969e-07 ~ 1,
      p >= 1.186969e-07 ~ 0
    )
  )

table(cg09117673_simu_268_2140_1000$sig)

# 0   1 
# 104 896 

write_csv(cg09117673_simu_268_2140_1000, here("ichd_followup", "results", "power", "cg09117673_tuk_268_2140_1000.csv"))

# nshd prep ----
info <- read_csv(here("ichd", "info", "info_NSHD99_ASCVD_107cases_670controls_cardinality.csv"))
meth <- fread(here("master_data", "datasets", "NSHD", "beta_99_enmix.csv"), data.table = FALSE)
rownames(meth) <- meth[, 1]
meth <- meth[, -1]
meth <- meth[, match(info$Basename, colnames(meth))]

factor_cols <- c("sex", "smok", "Sentrix_ID", "Sentrix_Position", "ASCVD", "Status")
info[factor_cols] <- lapply(info[factor_cols], factor)

# inspect original ewas results
load(here("master_data", "cpg_lists", "all_bl_sex_cg.RData"))
resu <- read_csv(here("ichd", "results", "Resu07_iASCVD_NSHD99_107cases_670controls_cardinality_LinFix_m1.csv"))
resu <- resu %>%
  anti_join(sex450, by = c("cpg" = "cg")) %>%
  anti_join(sexepic, by = c("cpg" = "cg")) %>%
  anti_join(bl450, by = c("cpg" = "cg")) %>%
  anti_join(blepic, by = c("cpg" = "cg"))
# 750119 cpgs
resu$fdr <- p.adjust(resu$p, method = "fdr")

# nshd cg09117673 ----
meth_cg09117673 <- meth %>% filter(rownames(meth) == "cg09117673")
meth_cg09117673 <- data.frame(t(meth_cg09117673))
info$cg09117673 <- meth_cg09117673$cg09117673
ggplot(info, aes(ASCVD, cg09117673)) +
  geom_boxplot()
info_case <- info %>% filter(ASCVD == "1")
info_control <- info %>% filter(ASCVD == "0")

mean(info_case$cg09117673, na.rm = TRUE) - mean(info_control$cg09117673, na.rm = TRUE)
# 0.001631358

## power analysis with seeds ----
set.seed(1701)
cg09117673_simu_107_670_1000 <- power_simulation(107, 670, 1000, 
                                                  case_df = info_case, control_df = info_control,
                                                  full_df = info,
                                                  barcode_colname = "Basename", cg_number = "cg09117673",
                                                  model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                  model_raw = cg09117673 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)

cg09117673_simu_107_670_1000 <- cg09117673_simu_107_670_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg09117673_simu_107_670_1000$sig)
#0   1 
#843 157 

write_csv(cg09117673_simu_107_670_1000, here("ichd_followup", "results", "power", "cg09117673_nshd_107_670_1000.csv"))


set.seed(1701)

cg09117673_simu_214_1340_1000 <- power_simulation(214, 1340, 1000, 
                                                 case_df = info_case, control_df = info_control,
                                                 full_df = info,
                                                 barcode_colname = "Basename", cg_number = "cg09117673",
                                                 model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                 model_raw = cg09117673 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)
cg09117673_simu_214_1340_1000 <- cg09117673_simu_214_1340_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg09117673_simu_214_1340_1000$sig)

# 0   1 
# 193 807 

write_csv(cg09117673_simu_214_1340_1000, here("ichd_followup", "results", "power", "cg09117673_nshd_214_1340_1000.csv"))


set.seed(1701)

cg09117673_simu_321_2010_1000 <- power_simulation(321, 2010, 1000, 
                                                  case_df = info_case, control_df = info_control,
                                                  full_df = info,
                                                  barcode_colname = "Basename", cg_number = "cg09117673",
                                                  model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                  model_raw = cg09117673 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)
cg09117673_simu_321_2010_1000 <- cg09117673_simu_321_2010_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg09117673_simu_321_2010_1000$sig)

# 0   1 
# 14 986 

write_csv(cg09117673_simu_321_2010_1000, here("ichd_followup", "results", "power", "cg09117673_nshd_321_2010_1000.csv"))



## explore incremental patterns

simu_cg09117673 <- function(case_n){
  df <- power_simulation(case_n, as.integer(case_n*670/107), 1000, 
                         case_df = info_case, control_df = info_control,
                         full_df = info,
                         barcode_colname = "Basename", cg_number = "cg09117673",
                         model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                         model_raw = cg09117673 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)
  
  df <- df %>%
    mutate(
      sig = case_when(
        p < 6.665609e-08 ~ 1,
        p >= 6.665609e-08 ~ 0
      )
    )
  
  pw <- c(case_n, as.integer(case_n*670/107), data.frame(table(df$sig))$Freq)
  
  pw
}


library(parallel) 
options(mc.cores = 6) 
resu <- do.call(rbind, mclapply(seq(50, 300, 10), simu_cg09117673))
resu <- data.frame(resu)
plot(resu$X1, resu$X4)
write_csv(resu, here("ichd_followup", "results", "power", "cg09117673_nshd_increment_noseed_1000.csv"))

set.seed(1701)
resu_seeded <- do.call(rbind, mclapply(seq(50, 300, 10), simu_cg09117673))
resu_seeded <- data.frame(resu_seeded)
plot(resu_seeded$X1, resu_seeded$X4)

# nshd cg26965354 ----
meth_cg26965354 <- meth %>% filter(rownames(meth) == "cg26965354")
meth_cg26965354 <- data.frame(t(meth_cg26965354))
info$cg26965354 <- meth_cg26965354$cg26965354
ggplot(info, aes(ASCVD, cg26965354)) +
  geom_boxplot()
info_case <- info %>% filter(ASCVD == "1")
info_control <- info %>% filter(ASCVD == "0")

mean(info_case$cg26965354, na.rm = TRUE) - mean(info_control$cg26965354, na.rm = TRUE)
# 0.002977642

set.seed(1701)

cg26965354_simu_107_670_1000 <- power_simulation(107, 670, 1000, 
                                                 case_df = info_case, control_df = info_control,
                                                 full_df = info,
                                                 barcode_colname = "Basename", cg_number = "cg26965354",
                                                 model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                 model_raw = cg26965354 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)

cg26965354_simu_107_670_1000 <- cg26965354_simu_107_670_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg26965354_simu_107_670_1000$sig)

# 0   1 
# 856 144 

write_csv(cg26965354_simu_107_670_1000, here("ichd_followup", "results", "power", "cg26965354_nshd_107_670_1000.csv"))


set.seed(1701)

cg26965354_simu_214_1340_1000 <- power_simulation(214, 1340, 1000, 
                                                 case_df = info_case, control_df = info_control,
                                                 full_df = info,
                                                 barcode_colname = "Basename", cg_number = "cg26965354",
                                                 model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                 model_raw = cg26965354 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)

cg26965354_simu_214_1340_1000 <- cg26965354_simu_214_1340_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg26965354_simu_214_1340_1000$sig)

#   0   1 
# 214 786

write_csv(cg26965354_simu_214_1340_1000, here("ichd_followup", "results", "power", "cg26965354_nshd_214_1340_1000.csv"))


set.seed(1701)

cg26965354_simu_321_2010_1000 <- power_simulation(321, 2010, 1000, 
                                                  case_df = info_case, control_df = info_control,
                                                  full_df = info,
                                                  barcode_colname = "Basename", cg_number = "cg26965354",
                                                  model_norm = norm ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position, 
                                                  model_raw = cg26965354 ~ ASCVD + BMI + CD8T + CD4T + NK + Bcell + Mono + Gran + sex + smok + Status + Sentrix_ID + Sentrix_Position)

cg26965354_simu_321_2010_1000 <- cg26965354_simu_321_2010_1000 %>%
  mutate(
    sig = case_when(
      p < 6.665609e-08 ~ 1,
      p >= 6.665609e-08 ~ 0
    )
  )
table(cg26965354_simu_321_2010_1000$sig)

# 0   1 
# 14 986 

write_csv(cg26965354_simu_321_2010_1000, here("ichd_followup", "results", "power", "cg26965354_nshd_321_2010_1000.csv"))
