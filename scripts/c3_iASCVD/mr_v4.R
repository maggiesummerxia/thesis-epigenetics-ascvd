library(here)
library(tidyverse)
library(readxl)
library(data.table)
library(MendelianRandomization)
library(ieugwasr)

# import ----
## import cpg sites needed ----

signal_gwas_rep <- c("cg04057093", "cg16325826", "cg20248056", "cg23337361", "cg10355862")
signal_ewas_rep <- c("cg22304262", "cg26467725", "cg00355799")

resu03 <- read_csv(here("ichd", "results", "v3", "meta_v3_p03_372cpg.csv"))
cpg <- tibble("cpg" = c(resu03$rs_number, signal_ewas_rep))
cpg <- cpg %>%
  mutate(comments =
           case_when(
             cpg %in% signal_gwas_rep ~ "gwas_rep_v3",
             cpg %in% signal_ewas_rep ~ "ewas_rep"
           ))

## import epic manifest ----

mani <- read_rds(here("master_data", "manifest", "array_manifest_epic.rds"))

## import chd gwas results ----

chd_gwas <- read_tsv(here("ichd_followup", "data", "phs001672_tcheandjieu2022_full_sum_stats", "phs001672.pha005196.txt"), skip = 20)

# import meqtl results

qtl_450 <- read.csv(here("ichd_followup", "results", "mr", "meqtl_for_top03_n_cg22304262_450k_106cpg.csv"))
qtl_epic <- read.csv(here("ichd_followup", "results", "mr", "meqtl_for_top03_epic_67cpg.csv"))

# MR 450k ----
# select clumped cpgs
qtl_450 <- filter(qtl_450, clumped == "TRUE")

# rename
mr_450_cpg <- qtl_450[c("cpg", "rsid", "beta_a1", "se", "pval", "allele1", "allele2")]
colnames(mr_450_cpg) <- c("cg", "rsid", "beta_cpg", "se_cpg", "p_cpg", "a1_cpg", "a2_cpg")

mr_450_gwas <- chd_gwas %>% filter(`SNP ID` %in% mr_450_cpg$rsid)
mr_450_gwas <- mr_450_gwas[c("SNP ID", "|&beta;|", "SE", "P-value", "Allele1", "Allele2")]
colnames(mr_450_gwas) <- c("rsid", "beta_chd", "se_chd", "p_chd", "a1_chd", "a2_chd")

mrdf_450 <- left_join(mr_450_cpg, mr_450_gwas, by = "rsid")

# remove cpg with only 1 snp
cpg_n_snp_450 <- data.frame(table(mrdf_450$cg))
cpg_n_snp_450 <- cpg_n_snp_450 %>% filter(Freq > 1) #26 cpgs
mrdf_450 <- mrdf_450 %>% filter(cg %in% cpg_n_snp_450$Var1)

# harmonies
mrdf_450 <- mrdf_450 %>%
  mutate(
    beta_chd_hm = case_when(
      a1_chd == a1_cpg & a2_chd == a2_cpg ~ beta_chd,
      a1_chd == a2_cpg & a2_chd == a1_cpg ~ -beta_chd,
    )
  )

mrdf_450_list <- mrdf_450 %>%
  group_by(cg) %>%
  group_split(.keep = TRUE) # list of 26

run_mr <- function(dat){
  mr_ob = mr_input(bx = dat$beta_cpg, 
                   bxse = dat$se_cpg, 
                   by = dat$beta_chd_hm, 
                   byse = dat$se_chd)
  ivw_ob <- mr_ivw(mr_ob)
  
  resu <- data.frame("cg" = dat[1,]$cg, "estimate" = ivw_ob$Estimate, "se" = ivw_ob$StdError, "p" = ivw_ob$Pvalue, "n_snp" = ivw_ob$SNPs)
}

mr_resu_450_list <- lapply(mrdf_450_list, run_mr)

mr_resu_450 <-bind_rows(mr_resu_450_list)