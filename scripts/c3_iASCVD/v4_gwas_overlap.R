library(here)
library(tidyverse)
library(tidylog)
library(readxl)
library(data.table)
library(topGO)
library(fs)
library(biomaRt)

gwas1 <- read_excel(here("master_data", "cpg_lists", "Tcheandjieu_2022_supp_tables.xlsx"), sheet = 3, skip = 1)
gwas1 <- gwas1[, c("rsID", "Chr", "Position", "region", "Effect", "StdErr", "P-value", "PMID", "PUBLISHED")]
colnames(gwas1) <- c("rsid_38", "snp_chr", "snp_pos", "gene", "snp_effect", "snp_se", "snp_p", "pmid", "publish_date")
gwas1_sig <- gwas1 %>% filter(snp_p < 5e-08) #36

gwas2 <- read_excel(here("master_data", "cpg_lists", "Tcheandjieu_2022_supp_tables.xlsx"), sheet = 4, skip = 1)
gwas2 <- gwas2[-34,]
gwas2 <- gwas2[, c("rsID", "chr", "pos", "Gene symbol", "beta", "SE", "P")]
colnames(gwas2) <- c("rsid_38", "snp_chr", "snp_pos", "gene", "snp_effect", "snp_se", "snp_p")
gwas2$pmid <- 35915156
gwas2$publish_date <- as.Date("2022-08-01", format = "%Y-%m-%d")

gwas_1n2 <- bind_rows(gwas1, gwas2)

gwas_lead <- gwas_1n2[, c("rsid_38", "pmid", "publish_date")]
gwas_lead$lead <- "Y"

tcheandjieu_all <- read_tsv(here("ichd_followup", "data", "phs001672_tcheandjieu2022_full_sum_stats", "phs001672.pha005196.txt"), skip = 20)
ol_all <- inner_join(tcheandjieu_all, gwas_lead, by = c("SNP ID" = "rsid_38"))
# 16 snps are not in the file. 180 are.

# select bins from 
lead_bins <- unique(ol_all$`Bin ID`) #170 bins. 10 bins have 2 lead SNPs.

tcheandjieu_170_bins <- tcheandjieu_all %>%
  filter(`Bin ID` %in% lead_bins)
# filter: removed 24,088 rows (37%), 40,459 rows remaining

# append lead snps
tcheandjieu_170_bins <- tcheandjieu_170_bins %>% left_join(gwas_lead, by = c("SNP ID" = "rsid_38"))
tcheandjieu_170_bins <- tcheandjieu_170_bins[, c("SNP ID", "P-value", "Chr ID", "Chr Position", 
                                                 "Allele1", "Allele2", "|&beta;|", "SE", "Sample size", 
                                                 "Bin ID", "pmid", "publish_date", "lead")]
colnames(tcheandjieu_170_bins) <- c("rsid", "gwas_p", "chr_38", "pos_38", 
                                    "gwas_allele1", "gwas_allele2", "gwas_beta", "gwas_se", "gwas_n", 
                                    "gwas_bin", "replicate_from_pmid", "publish_date", "gwas_lead")
# fill bins with pmid
bin_pmid <- data.frame(table(tcheandjieu_170_bins$gwas_bin, tcheandjieu_170_bins$replicate_from_pmid)) %>% 
  filter(Freq > 0) %>%
  arrange(desc(Freq))
# filter: removed 2,882 rows (94%), 178 rows remaining
# with only 170 unique snps, some are led by more then one study

tcheandjieu_170_bins <- tcheandjieu_170_bins %>%
  arrange(gwas_bin, gwas_lead, publish_date)

tcheandjieu_170_bins <- tcheandjieu_170_bins %>%
  group_by(gwas_bin) %>%
  mutate(pmid_first_publish = replicate_from_pmid[1]) #40455 rsid

# write_csv(tcheandjieu_170_bins, here("ichd_followup", "results", "gwas_overlap", "tcheandjieu_sig_170_bins.csv"))


# cpg overlap with snps
# overlap ewas cpgs with snps
qtl <- read_csv(here("ichd_followup", "results", "mr", "meqtl_for_top03_450k_epic_172cpg.csv"))
qtl_tcheandjieu_ol <- inner_join(qtl, tcheandjieu_170_bins, by = c("rsid" = "rsid"))
# 1,526 rows, 9 bins

# select only cis

godmc_full <- fread(here("master_data", "godmc", "assoc_meta_all.csv.gz"), data.table = FALSE)
godmc_snp <- fread(here("master_data", "godmc", "snps.csv.gz"), data.table = FALSE)
godmc_full_cistrans <- godmc_full[, c("cpg", "snp", "cistrans")]
godmc_snp <- godmc_snp[, c("name", "rsid")]
godmc_full_cistrans <- godmc_full_cistrans %>% inner_join(godmc_snp, by = c("snp" = "name"))
godmc_full_cistrans <- godmc_full_cistrans %>%
  dplyr::select(cpg, rsid, cistrans)

qtl_tcheandjieu_ol <- qtl_tcheandjieu_ol %>%
  left_join(godmc_full_cistrans, by = c("cpg"="cpg", "rsid"="rsid")) #1526
qtl_tcheandjieu_ol$cistrans[qtl_tcheandjieu_ol$meqtl_source == "sergio"] <- "TRUE"

# select only cis meQTL
qtl_tcheandjieu_ol <- qtl_tcheandjieu_ol %>% filter(cistrans == "TRUE")
# filter: removed 67 rows (4%), 1,459 rows remaining


# annotate rsid with genes
anno_xplore <- read_tsv(here("ichd_followup", "data", "v4_gwas_followup", "RESULTS_95752", "snp_annotation.txt"))
anno_xplore_gene <- anno_xplore[, c("ID", "snp_conseq", "snp_conseq_gene", "eqtl", "sqtl", "positional_mapping", "source_finalGenes", "geneList")]

qtl_tcheandjieu_ol <- qtl_tcheandjieu_ol %>% left_join(anno_xplore_gene, by = c("rsid" = "ID"))

cpg_n_bin <- data.frame(table(qtl_tcheandjieu_ol$cpg, qtl_tcheandjieu_ol$gwas_bin)) %>%
  filter(Freq > 0) # 10

qtl_tcheandjieu_ol <- qtl_tcheandjieu_ol %>% arrange(gwas_lead, gwas_p)

qtl_tcheandjieu_ol_top <- qtl_tcheandjieu_ol %>%
  group_by(cpg, gwas_bin) %>%
  filter(row_number()==1)

# group_by: 2 grouping variables (cpg, gwas_bin)
# filter (grouped): removed 1,451 rows (99%), 8 rows remaining

write_csv(qtl_tcheandjieu_ol, here("ichd_followup", "results", "gwas_overlap", "tcheandjieu_meqtl_overlap_8_bins.csv"))

write_csv(qtl_tcheandjieu_ol_top, here("ichd_followup", "results", "gwas_overlap", "gwas_overlap_table_for_paper.csv"))

# TODO: convert meQTL 37 to 38


# write_csv(gwas_1n2, here("ichd_followup", "data", "v4_gwas_followup", "35915156_Tcheandjieu_2022", "gwas_all_sig.csv"))

anno_xplore_gwas_sig <- read_tsv(here("ichd_followup", "data", "v4_gwas_followup", "RESULTS_35065", "snp_annotation.txt"))
anno_xplore_gwas_sig <- anno_xplore_gwas_sig[, c("ID", "snp_conseq", "snp_conseq_gene", "eqtl", "sqtl", "positional_mapping", "source_finalGenes", "geneList")]


# remove the pmids not involved in the big comparison


# TODO: annotate with biomart
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

gene_table <- getBM(attributes=c(
  "refsnp_id", "ensembl_gene_stable_id", "associated_gene", "ensembl_gene_name"),
                    filters = "snp_filter", values = gwas_1n2$rsid_38, mart = ensembl, useCache = FALSE)

gene_table <- gene_table %>% arrange(associated_gene)

gene_table_narm <- data.frame(gene_table) %>%
  filter(!associated_gene == "") %>%
  filter(!associated_gene == "intergenic")
# filter: removed 118 rows (27%), 327 rows remaining 
# filter: removed 37 rows (11%), 290 rows remaining

qtl_tcheandjieu_ol_top <- read_csv(here("ichd_followup", "results", "gwas_overlap", "gwas_overlap_table_for_paper.csv"))

qtl_ol_top_anno_info <- getBM(attributes=c(
  "refsnp_id", "ensembl_gene_stable_id", "associated_gene", "ensembl_gene_name"),
  filters = "snp_filter", values = qtl_tcheandjieu_ol_top$rsid, mart = ensembl, useCache = FALSE)



# TODO: redo gene overlap with new gwas gene list

# remove SNPs that are not involved in the mQTL overlap analysis
gwas_snp_final <- tcheandjieu_170_bins %>% drop_na(gwas_lead) #180 lead snps, 170 bins
gene_table_narm <- inner_join(gene_table_narm, gwas_snp_final["rsid"], by = c("refsnp_id" = "rsid")) 
# 147 snps are annotated to genes

# make gene list for 180 gwas signatures
gene_vector <- as.vector(gene_table_narm$associated_gene)
gene_concentrate <- paste(gene_vector, collapse=",")
gene_vector_tidy <- strsplit(gene_concentrate, ",")[[1]]

gwas_gene_list <- data.frame(X1 = gene_vector_tidy)
gwas_gene_list <- unique(gwas_gene_list) #325 genes

ewas_resu <- read_csv(here("ichd", "results", "v3", "meta_resu_v3_p03_372cpg.csv"))
ewas_resu <- ewas_resu %>% drop_na("UCSC_RefGene_Name")
# drop_na: removed 117 rows (31%), 255 rows remaining

ewas_gene_list <- read_tsv(here("ichd", "results", "v3", "meta_v3_p03_372cpg_gene_list.txt"), col_names = FALSE) 
ewas_gene_list <- unique(ewas_gene_list) #279 genes

# overlap
gene_ol <- inner_join(gwas_gene_list, ewas_gene_list, by = "X1")
# rows total 5

gene_ol <- gene_ol$X1

extract_cpg_contains_certain_gene <- function(gene, dat){
  dat[dat$UCSC_RefGene_Name %like% gene, ]
}

gwas_ewas_ol <- lapply(gene_ol, extract_cpg_contains_certain_gene, dat = ewas_resu) %>% bind_rows()
gwas_ewas_ol <- gwas_ewas_ol %>% arrange(`p-value`)


# check SNP build GRCh38/hg38 or GRCh37/hg19
## godmc
 # hg19/build37
## sergio's meqtl zz
 # GRCh37/hg19
## tcheandjieu
 # Human genome build:	38
 # dbSNP build:	141
## van der harst
 # 29212778-GCST005195-EFO_0000378.h.tsv.gz is 38
## epic menifest
 # 37
