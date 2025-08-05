library(here)
library(tidyverse)
library(QCEWAS)

# change col names for GS results
age <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet1_Menopause-EWAS_20220228/MENOPAUSE_GENSCOT1_MODELA_DLM_21022022.csv")
colnames(age) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(age, here("menopause", "results", "age", "age_GS1.txt"))

age2 <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet2_Menopause-EWAS_20220228/MENOPAUSE_GENSCOT2_MODELA_DLM_21022022.csv")
colnames(age2) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(age2, here("menopause", "results", "age", "age_GS2.txt"))

colnames(duration) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(duration, here("menopause", "results", "duration", "duration_GS1.txt"))

duration2 <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet2_Menopause-EWAS_20220228/DURATION_GENSCOT2_MODELA_DLM_21022022.csv")
colnames(duration2) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(duration2, here("menopause", "results", "duration", "duration_GS2.txt"))

status <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet1_Menopause-EWAS_20220228/MENOPAUSALSTATUS_GENSCOT1_MODELA_DLM_21022022.csv")
colnames(status) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(status, here("menopause", "results", "status", "status_GS1.txt"))

status2 <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet2_Menopause-EWAS_20220228/MENOPAUSALSTATUS_GENSCOT2_MODELA_DLM_21022022.csv")
colnames(status2) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(status2, here("menopause", "results", "status", "status_GS2.txt"))

hrt <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet1_Menopause-EWAS_20220228/HRTUSE_GENSCOT1_MODELA_DLM_21022022.csv")
colnames(hrt) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(hrt, here("menopause", "results", "hrt", "hrt_GS1.txt"))

hrt2 <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet2_Menopause-EWAS_20220228/HRTUSE_GENSCOT2_MODELA_DLM_21022022.csv")
colnames(hrt2) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(hrt2, here("menopause", "results", "hrt", "hrt_GS2.txt"))

timing <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet1_Menopause-EWAS_20220228/EARLYLATEMENOPAUSE_GENSCOT1_MODELA_DLM_21022022.csv")
colnames(timing) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(timing, here("menopause", "results", "timing", "timing_GS1.txt"))

timing2 <- read_csv("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Projects/Menopause/GS/GenScotSet2_Menopause-EWAS_20220228/EARLYLATEMENOPAUSE_GENSCOT2_MODELA_DLM_21022022.csv")
colnames(timing2) <- c("MARKER", "CHR", "POS", "BETA", "SE", "TSTATISTIC", "PVAL", "FDR", "N")
write_tsv(timing2, here("menopause", "results", "timing", "timing_GS2.txt"))

# preparation
translation <- data_frame("ori" = c("PROBEID", "P_VAL"), "new" = c("MARKER", "PVAL"))
# translation_meta <- data_frame("ori" = c("PROBEID", "BETA", "SE", "N","P_VAL"), "new" = c("rs_number", "beta", "se", "n_samples","`p-value`"))
load(here("master_data", "cpg_lists", "all_bl_sex_cg.RData"))

# menopause age
setwd(here("menopause", "results", "age"))

age_results_uk <- c("age_TWINSUK.txt", "age_NSHD99.txt", "age_NSHD09ken.txt", 
                    "age_NSHD09bell.txt", "age_NCDS.txt", "age_GS1.txt",
                    "age_GS2.txt", "age_EPICNORFOLK.txt", "age_ALSPAC.txt")
QC_results_age <- EWAS_series(EWAS_files = age_results_uk,
                          output_files = c("TwinsUK", "NSHD1", "NSHD2", "NSHD3", "NCDS", "GS1", "GS2", "EPICNORFOLK", "ALSPAC"),
                          save_final_dataset = FALSE,
                          header_translations = translation)

meta_age <- read_table(here("menopause", "results", "age", "gwama_uk8", "gwama_uk8.out"))
meta_age <- meta_age %>% 
  anti_join(bl450, by = c("rs_number" = "cg")) %>%
  anti_join(blepic, by = c("rs_number" = "cg")) %>%
  anti_join(sex450, by = c("rs_number" = "cg")) %>%
  anti_join(sexepic, by = c("rs_number" = "cg"))
meta_age$fdr <- p.adjust(meta_age$`p-value`, method = "fdr")
# write_csv(meta_age, here("menopause", "results", "age", "gwama_uk8", "gwama_uk8_age_blrm_sexrm.csv"))

# reproductive duration
setwd(here("menopause", "results", "duration"))

duration_results_uk <- c("duration_TWINSUK.txt", "duration_NSHD99.txt", "duration_NSHD09ken.txt", 
                    "duration_NSHD09bell.txt", "duration_NCDS.txt", "duration_GS1.txt",
                    "duration_GS2.txt", "duration_EPICNORFOLK.txt", "duration_ALSPAC.txt")
QC_results_duration <- EWAS_series(EWAS_files = duration_results_uk,
                          output_files = c("TwinsUK", "NSHD1", "NSHD2", "NSHD3", "NCDS", "GS1", "GS2", "EPICNORFOLK", "ALSPAC"),
                          save_final_dataset = FALSE,
                          header_translations = translation)

meta_duration <- read_table(here("menopause", "results", "duration", "gwama_uk8", "gwama_uk8.out"))
meta_duration <- meta_duration %>% 
  anti_join(bl450, by = c("rs_number" = "cg")) %>%
  anti_join(blepic, by = c("rs_number" = "cg")) %>%
  anti_join(sex450, by = c("rs_number" = "cg")) %>%
  anti_join(sexepic, by = c("rs_number" = "cg"))
meta_duration$fdr <- p.adjust(meta_duration$`p-value`, method = "fdr")

# write_csv(meta_duration, here("menopause", "results", "duration", "gwama_uk8", "gwama_uk8_duration_blrm_sexrm.csv"))

setwd(here("menopause", "results", "duration", "gwama_uk8"))
meta_duration <- meta_duration %>% rename("P_VAL" = "p-value")
EWAS_plots(meta_duration,
           plot_QQ = TRUE,
           plot_Man = FALSE,
           plot_cutoff_p = 0.05,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = "gwama_uk8_duration")

# timing
setwd(here("menopause", "results", "timing"))

timing_results_uk <- c("timing_TWINSUK.txt", "timing_NSHD09ken.txt", 
                         "timing_NSHD09bell.txt", "timing_GS1.txt",
                         "timing_GS2.txt", "timing_EPICNORFOLK.txt")
QC_results_timing <- EWAS_series(EWAS_files = timing_results_uk,
                                   output_files = c("TwinsUK", "NSHD2", "NSHD3", "GS1", "GS2", "EPICNORFOLK"),
                                   save_final_dataset = FALSE,
                                   header_translations = translation)

meta_timing <- read_table(here("menopause", "results", "timing", "gwama_uk6", "gwama_uk6.out"))
meta_timing <- meta_timing %>% 
  anti_join(bl450, by = c("rs_number" = "cg")) %>%
  anti_join(blepic, by = c("rs_number" = "cg")) %>%
  anti_join(sex450, by = c("rs_number" = "cg")) %>%
  anti_join(sexepic, by = c("rs_number" = "cg"))
meta_timing$fdr <- p.adjust(meta_timing$`p-value`, method = "fdr")

# write_csv(meta_timing, here("menopause", "results", "timing", "gwama_uk6", "gwama_uk6_timing_blrm_sexrm.csv"))

setwd(here("menopause", "results", "timing", "gwama_uk6"))
meta_timing <- meta_timing %>% rename("P_VAL" = "p-value")
EWAS_plots(meta_timing,
           plot_QQ = TRUE,
           plot_Man = FALSE,
           plot_cutoff_p = 0.05,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = "gwama_uk6_timing")

# status
setwd(here("menopause", "results", "status"))

status_results_uk <- c("status_TWINSUK.txt", "status_NSHD99.txt", "status_NCDS.txt", "status_GS1.txt",
                         "status_GS2.txt", "status_EPICNORFOLK.txt", "status_ALSPAC.txt")
QC_results_status <- EWAS_series(EWAS_files = status_results_uk,
                                   output_files = c("TwinsUK", "NSHD1", "NCDS", "GS1", "GS2", "EPICNORFOLK", "ALSPAC"),
                                   save_final_dataset = FALSE,
                                   header_translations = translation)

meta_status <- read_table(here("menopause", "results", "status", "gwama_uk6", "gwama_uk6.out"))
meta_status <- meta_status %>% 
  anti_join(bl450, by = c("rs_number" = "cg")) %>%
  anti_join(blepic, by = c("rs_number" = "cg")) %>%
  anti_join(sex450, by = c("rs_number" = "cg")) %>%
  anti_join(sexepic, by = c("rs_number" = "cg"))
meta_status$fdr <- p.adjust(meta_status$`p-value`, method = "fdr")

# write_csv(meta_status, here("menopause", "results", "status", "gwama_uk6", "gwama_uk6_status_blrm_sexrm.csv"))

setwd(here("menopause", "results", "status", "gwama_uk6"))
meta_status <- meta_status %>% rename("P_VAL" = "p-value")
EWAS_plots(meta_status,
           plot_QQ = TRUE,
           plot_Man = FALSE,
           plot_cutoff_p = 0.05,
           plot_QQ_bands = FALSE,
           high_quality_plots = FALSE,
           save_name = "gwama_uk6_status")

# hrt
setwd(here("menopause", "results", "hrt"))

hrt_results_uk <- c("hrt_TWINSUK.txt", "hrt_NSHD99.txt", "hrt_GS1.txt",
                         "hrt_GS2.txt", "hrt_EPICNORFOLK.txt")
QC_results_hrt <- EWAS_series(EWAS_files = hrt_results_uk,
                                   output_files = c("TwinsUK", "NSHD1", "GS1", "GS2", "EPICNORFOLK"),
                                   save_final_dataset = FALSE,
                                   header_translations = translation)

meta_hrt <- read_table(here("menopause", "results", "hrt", "gwama_uk5", "gwama_uk5.out"))
meta_hrt <- meta_hrt %>% 
  anti_join(bl450, by = c("rs_number" = "cg")) %>%
  anti_join(blepic, by = c("rs_number" = "cg")) %>%
  anti_join(sex450, by = c("rs_number" = "cg")) %>%
  anti_join(sexepic, by = c("rs_number" = "cg"))
meta_hrt$fdr <- p.adjust(meta_hrt$`p-value`, method = "fdr")

# write_csv(meta_hrt, here("menopause", "results", "hrt", "gwama_uk5", "gwama_uk5_hrt_blrm_sexrm.csv"))

# replication of signals in US and EU
duration_us <- read_table("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/A2_Duration/gwama_us2/gwama_us2.in.out")
duration_eu <- read_table("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/A2_Duration/gwama_eu2/gwama_eu2.out")

timing_us <- read_table("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/A3_EarlyLate/gwama_us2/gwama_us2.out")
timing_eu <- read_table("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/A3_EarlyLate/gwama_eu2/gwama_eu2.out")

status_eu <- read_table("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/A4_Status/gwama_eu2/gwama_eu2.out")

hrt_twin <- read_csv("/Users/maggiexia/OneDrive/Projects/Menopause/Useful_results/discMZ/gwama_mz_trimmed_random_new/fdr5_mz_new_random.csv")

# forest plots