# iASCVD paper v4 make b vitamins plots
library(tidyverse)
library(here)
library(data.table)

meth <- fread(here("master_data", "datasets", "TwinsUK", "450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv"), data.table = FALSE)

cg22304262 <- meth %>% filter(cg == "cg22304262")
cg22304262 <- data.frame(t(cg22304262))
cg22304262$Barcode <- rownames(cg22304262)
colnames(cg22304262) <- c("cg22304262", "Barcode")
cg22304262 <- cg22304262[-1, ]

info <- read_csv(here("master_data", "info", "twinsuk", "info_450k_date_gathered.csv"))
info <- info[, c("KCLid", "Barcode")]

cg22304262 <- left_join(cg22304262, info, by = "Barcode")

bv <- read_csv(here("ichd_followup", "data", "diet", "b12_biotin_pantothenate_487ids.csv")) #487

bv <-left_join(bv, cg22304262, by = "KCLid")
bv$cg22304262 <- as.numeric(bv$cg22304262)

bv %>%
  ggplot(aes(x = cg22304262, y = adj_vitamin_b12_ug)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x, se = F) +
  xlab("cg22304262 methylation level") +
  ylab("Vitamin B12 level") +
  theme_bw()

bv %>%
  ggplot(aes(x = cg22304262, y = adj_pantothenate_mg)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x, se = F) +
  xlab("cg22304262 methylation level") +
  ylab("Pantothenate level") +
  theme_bw()

bv %>%
  ggplot(aes(x = cg22304262, y = adj_biotin_ug)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x, se = F) +
  xlab("cg22304262 methylation level") +
  ylab("Biotin level") +
  theme_bw()

# manually export 300*300