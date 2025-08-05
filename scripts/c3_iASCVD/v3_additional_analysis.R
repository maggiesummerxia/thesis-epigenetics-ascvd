library(here)
library(tidyverse)
library(readxl)
library(tidylog)

# prep ----

resu <- read_csv(here("ichd", "results", "v3", "meta_v3_p03_372cpg.csv"))
resu <- resu %>%
  arrange("p-value") %>%
  mutate(rank = 1:372)
resu_short <- resu[c("rs_number", "array_label", "rank", "beta", "se", "p-value", "fdr", 
                     "n_samples", "effects", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
colnames(resu_short) <- c("cpg", "cpg_array_label", "cpg_rank", "cpg_beta", "cpg_se", "cpg_p-value", "cpg_fdr", 
                          "cpg_n_samples", "cpg_effects", "cpg_chr", "cpg_pos", "cpg_gene", "cpg_group")

mqtl_450 <- read_csv(here("ichd_followup", "results", "mr", "meqtl_for_top03_450k_105cpg.csv"))
mqtl_450_list <- mqtl_450[c("cpg", "rsid")] #105 cpg, 25551 rows
mqtl_450_list$meqtl_source <- "godmc"
# warning: this includes cg22304262, so among the top 372, only 105 has meQTL

mqtl_epic <- read_csv(here("ichd_followup", "results", "mr", "meqtl_for_top03_epic_67cpg.csv"))
mqtl_epic_list <- mqtl_epic[c("CpG", "rsid")] # 67 cpg, 8509 rows
colnames(mqtl_epic_list) <- c("cpg", "rsid")
mqtl_epic_list$meqtl_source <- "sergio"

mqtl_epic_clumped <- read_rds(here("ichd_followup", "results", "mr", "clump", "mqtl_epic_clump_2000_01_1e5.rds"))

# some of them become index so can't combine the clumped + unclumped

mqtl <- bind_rows(mqtl_450_list, mqtl_epic_list) # 172 unique cpgs

# annotate mqtl with resu_short
mqtl <- left_join(mqtl, resu_short, by = "cpg")

signal_gwas_rep <- c("cg04057093", "cg16325826", "cg20248056", "cg23337361", "cg10355862")
mqtl <- mqtl %>%
  mutate(comments =
           case_when(
             cpg %in% signal_gwas_rep ~ "gwas_rep"
           ))

# write_csv(mqtl, here("ichd_followup", "results", "mr", "meqtl_for_top03_450k_epic_172cpg.csv"))

gwas1 <- read_excel(here("master_data", "cpg_lists", "Tcheandjieu_2022_supp_tables.xlsx"), sheet = 3, skip = 1)
gwas1 <- gwas1[c("rsID", "Chr", "Position", "region")]
colnames(gwas1) <- c("rsid", "snp_chr", "snp_pos", "snp_gene")
gwas2 <- read_excel(here("master_data", "cpg_lists", "Tcheandjieu_2022_supp_tables.xlsx"), sheet = 4, skip = 1)
gwas2 <- gwas2[c("rsID", "chr", "pos", "Gene symbol")]
colnames(gwas2) <- c("rsid", "snp_chr", "snp_pos", "snp_gene")

gwas_list_white <- drop_na(bind_rows(gwas1, gwas2)) #196 rsid

# overlap between meqtl SNPs and the big GWAS hits ----

overlap_mqtl_gwas <- inner_join(mqtl, gwas_list_white, by = "rsid")
# write_csv(overlap_mqtl_gwas, here("ichd_followup", "results", "mr", "meqtl_gwas_overlap_for_top.csv"))

# what about the GWAS replication ones
mqtl_rep <- mqtl_450 %>% filter(cpg %in% signal_gwas_rep) # 319 rows
table(mqtl_rep$cpg, mqtl_rep$clumped) 
# 1 for cg16325826 and 318 for cg20248056
# 2 for cg20248056 after clumping

# overlap between historical EWAS hits and the big GWAS hits ----

gwas_gene <- c(gwas1$snp_gene, gwas2$snp_gene)
gwas_gene <- na.omit(gwas_gene)

ewas_his <- read_csv(here("master_data", "cpg_lists", "Candidate_gene_list_all.csv"))
ewas_his <- ewas_his %>% filter(Group == "iCHD")
ewas_gene <- na.omit(c(ewas_his$Gene1, ewas_his$Gene2))

intersect(gwas_gene, ewas_gene)

# annotate MR results ----
mr_resu <- read_csv(here("ichd_followup", "results", "mr", "v3_resu_mr_450_epic_m1.csv"))
mr_resu <- left_join(mr_resu, resu_short, by = c("cg" = "cpg")) # no info for cg22304262
mr_resu <- mr_resu %>% 
  drop_na(se) %>% #53 rows remaining 
 arrange(p)

mr_resu[c("estimate", "se", "cpg_beta", "cpg_se", "cpg_fdr" )] <- round(mr_resu[c("estimate", "se", "cpg_beta", "cpg_se", "cpg_fdr" )], 3)


write_csv(mr_resu, here("ichd_followup", "results", "mr", "v3_resu_mr_450_epic_m1_anno_narm.csv"))
