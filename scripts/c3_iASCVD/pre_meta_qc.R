library(here)
library(tidyverse)
library(QCEWAS)

# post-ewas qc
# read in results
twinsuk <- read_tsv(here("ichd", "results", "Tresu01_iASCVD_LinFix_m1.txt"))
nshd <- read_tsv(here("ichd", "results", "nshd_ichd_resu_for_meta.txt"))
ncds <- read_tsv(here("ichd", "results", "ncds_ichd_resu_for_meta.txt"))

twinsuk <- twinsuk[, 1:5]
colnames(twinsuk) <- c("PROBEID", "BETA", "SE", "P_VAL", "N")
colnames(nshd) <- c("PROBEID", "BETA", "SE", "P_VAL", "N")
colnames(ncds) <- c("PROBEID", "BETA", "SE", "P_VAL", "N")

# removal of black lists and sex chr
bl450 <- read_csv(here("master_data", "cpg_lists", "450k_exclusion_probes_ids.csv"))
blepic <- read_csv(here("master_data", "cpg_lists", "EPIC_exclusion_probes_2020_update_v1_no_Pvars.csv"))
sex450 <- read_csv(here("master_data", "cpg_lists", "450k_sex_chr.csv"))
sexepic <- read_csv(here("master_data", "cpg_lists", "EPIC_sex_chr_cg.csv"))

sexepic <- sexepic[, 1]
colnames(bl450) <- "cg"
colnames(sex450) <- "cg"

twinsuk <- twinsuk %>% 
  anti_join(bl450, by = c("PROBEID" = "cg")) %>%
  anti_join(sex450, by = c("PROBEID" = "cg"))
# no change

nshd <- nshd %>% 
  anti_join(blepic, by = c("PROBEID" = "cg")) %>%
  anti_join(sexepic, by = c("PROBEID" = "cg"))
# no change

ncds <- ncds %>% 
  anti_join(blepic, by = c("PROBEID" = "cg")) %>%
  anti_join(sexepic, by = c("PROBEID" = "cg"))
# no change

#write_tsv(twinsuk, here("ichd", "results", "twinsuk_ichd_resu_for_qcewas"))
#write_tsv(nshd, here("ichd", "results", "nshd_ichd_resu_for_qcewas"))
#write_tsv(ncds, here("ichd", "results", "ncds_ichd_resu_for_qcewas"))

# multiple QCEWAS
setwd(here("ichd", "results"))
ichd_results <- c("twinsuk_ichd_resu_for_qcewas", "nshd_ichd_resu_for_qcewas", "ncds_ichd_resu_for_qcewas")
# QC_results <- EWAS_series(EWAS_files = ichd_results,
#                           output_files = c("TwinsUK", "NSHD", "NCDS"),
#                           save_final_dataset = FALSE)

# QCEWAS for meta analysis result
meta <- read_table(here("ichd", "results", "meta_analysis", "gwama_3cohorts", "gwama_3cohorts_random.out"))
meta$FDR <- p.adjust(meta$`p-value`, method = "fdr")
meta <- meta %>%
  rename("PROBEID" = "rs_number") %>%
  rename("BETA" = "beta") %>%
  rename("SE" = "se") %>%
  rename("P_VAL" = "p-value") %>%
  rename("N" = "n_samples")

setwd(here("ichd", "results", "meta_analysis", "gwama_3cohorts"))
# QC_results_meta <- EWAS_QC(data = meta,
#                            outputname = "Meta-analysis",
#                            save_final_dataset = FALSE)

# signal replication
meta_slim <- meta[, c("PROBEID", "BETA", "SE", "P_VAL", "FDR", 
                      "q_p-value", "i2", "n_studies", "N", "effects")]
list_ichd <- read_csv(here("master_data", "cpg_lists", "Candidate_gene_list_all.csv"))
list_lipid <- read_csv(here("master_data", "cpg_lists", "all_lipid_signals.csv"))
list_lip_mg <- read_csv(here("master_data", "cpg_lists", "maggie_lipid_bon_sig_random_signals.csv"))

list_ichd <- list_ichd %>% left_join(meta_slim, by = c("cpg" = "PROBEID"))
list_lipid <- list_lipid %>% left_join(meta_slim, by = c("MarkerName" = "PROBEID"))
list_lip_mg <- list_lip_mg %>% left_join(meta_slim, by = c("rs_number" = "PROBEID"))

list_ichd <- drop_na(list_ichd, BETA) #2018
list_lipid <- drop_na(list_lipid, BETA) #1340
list_lip_mg <- drop_na(list_lip_mg, BETA) #20

list_ichd_sig <- filter(list_ichd, P_VAL < 0.05)
list_lipid_sig <- filter(list_lipid, P_VAL < 0.05)
list_lip_mg_sig <- filter(list_lip_mg, P_VAL < 0.05)
