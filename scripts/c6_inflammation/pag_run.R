library(here)

meth_name <- here("master_data", "datasets", "NCDS1958", "Beta_enmix_TwinsUKnDimension_blrm_643barcode_777103cg.csv")
info_name <- here("pag", "data", "pag_paa_tuk_predict_combined_160id.csv")

f1 <- y ~ norm_all_100000011_phenylacetate + 
  Age + BMI + date_diff + 
  Smoking + Batch + STUDY_NAME + 
  CD8T + CD4T + NK + Bcell + Mono + Gran + 
  (1|Sentrix_ID) + (1|Sentrix_Position) + (1|KCLfam) + (1|Zygosity)

f2 <- y ~  Age + BMI + date_diff + 
  Smoking + Batch + STUDY_NAME + 
  CD8T + CD4T + NK + Bcell + Mono + Gran + 
  (1|Sentrix_ID) + (1|Sentrix_Position) + (1|KCLfam) + (1|Zygosity)

resu_name <- here("pag", "results", "resu_paa_value.csv")

source(here("pag", "scripts", "ewas_combined_epic_master.R"))
