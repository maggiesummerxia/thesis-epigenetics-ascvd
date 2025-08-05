library(fs) # for file system operations
library(here) # for hiding path names
library(tidyverse) # general data wrangling and visualisation
library(tidylog) # feedback about dplyr and tidyr operations
library(meta)
library(QCEWAS)

##
menosubset <- function(df, cpg){
  df <- df %>% 
    filter(MARKER %in% cpg) %>%
    select(STUDY, MARKER, BETA, SE, N) %>%
    mutate(SD = SE*(N^0.5))
  return(df)}

##
menoforest <- function(cpg_group){
  metaob <- metamean(n = cpg_group$N, 
                     mean = cpg_group$BETA, 
                     sd = cpg_group$SD, 
                     studlab = cpg_group$STUDY,
                     fixed = FALSE)
  forest(metaob)
}

##
change_class_chr <- function(df) {
  within(df, CHR <- as.character(CHR))
}


### duration
duration_data_path <- here("menopause", "results", "duration")

duration_results <- dir_ls(duration_data_path, regexp = "\\.txt$")

duration_data_list <- duration_results %>%
  map(~ read_tsv(
    .,
    col_names = TRUE,
    trim_ws = TRUE,
    name_repair = "unique"
  ))

# duration_data_names_raw <- names(duration_data_list)

duration_data_names <- duration_data_list %>%
  names() %>%
  str_remove_all(., paste0(duration_data_path)) %>%
  str_remove_all(., "/") %>%
  str_remove(., ".txt")

names(duration_data_list) <- duration_data_names

###
duration_data_list <- duration_data_list[c("duration_ALSPAC", "duration_TWINSUK", "duration_NSHD1", "duration_NSHD2", "duration_NSHD3",
                                           "duration_NCDS", "duration_GS1", "duration_GS2", "duration_EPICNORFOLK", "duration_FHS", 
                                           "duration_WHI", "duration_KORA", "duration_ROTTERDAM")]

duration_data_list <- lapply(duration_data_list, change_class_chr)

duration_data_all <- bind_rows(duration_data_list, .id = "STUDY")
duration_data_all$STUDY<-gsub("duration_", "", as.character(duration_data_all$STUDY))
ggplot(duration_data_all, aes(x = STUDY, y = SE)) + geom_boxplot()
ggplot(duration_data_all, aes(x = STUDY, y = BETA)) + geom_boxplot()

  # subset duration_data_all to produce boxplot for thesis






###

duration_cpg <- c("cg25946459", "cg07399807", "cg04794268", "cg13709442")

duration_data_all <- menosubset(duration_data_all, duration_cpg)

cg25946459 <- filter(duration_data_all, MARKER == "cg25946459")
cg25946459_uk <- filter(cg25946459, STUDY %in% c("TWINSUK", "NSHD1", "NSHD2", "NSHD3", "NCDS", 
                                              "GS1", "GS2", "EPICNORFOLK"))
cg25946459_us <- filter(cg25946459, STUDY %in% c("FHS", "WHI"))

cg25946459_eu <- filter(cg25946459, STUDY %in% c("KORA", "ROTTERDAM"))
png(here("menopause", "figures", "forest_duration_cg25946459_all.png"), width = 650, height = 350)
menoforest(cg25946459)
dev.off() 

png(here("menopause", "figures", "forest_duration_cg25946459_uk.png"), width = 650, height = 350)
menoforest(cg25946459_uk)
dev.off() 

menoforest(cg25946459_us)
menoforest(cg25946459_eu)

cg07399807 <- filter(duration_data_all, MARKER == "cg07399807")
cg07399807_uk <- filter(cg07399807, STUDY %in% c("TWINSUK", "NSHD1", "NSHD2", "NSHD3", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg07399807_us <- filter(cg07399807, STUDY %in% c("FHS", "WHI"))
cg07399807_eu <- filter(cg07399807, STUDY %in% c("KORA", "ROTTERDAM"))

png(here("menopause", "figures", "forest_duration_cg07399807_uk.png"), width = 650, height = 350)
menoforest(cg07399807_uk)
dev.off() 

cg04794268 <- filter(duration_data_all, MARKER == "cg04794268")
cg04794268_uk <- filter(cg04794268, STUDY %in% c("TWINSUK", "NSHD1", "NSHD2", "NSHD3", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg04794268_us <- filter(cg04794268, STUDY %in% c("FHS", "WHI"))
cg04794268_eu <- filter(cg04794268, STUDY %in% c("KORA", "ROTTERDAM"))

png(here("menopause", "figures", "forest_duration_cg04794268_all.png"), width = 650, height = 350)
menoforest(cg04794268)
dev.off() 

png(here("menopause", "figures", "forest_duration_cg04794268_uk.png"), width = 650, height = 350)
menoforest(cg04794268_uk)
dev.off() 

png(here("menopause", "figures", "forest_duration_cg04794268_us.png"), width = 650, height = 350)
menoforest(cg04794268_us)
dev.off() 

png(here("menopause", "figures", "forest_duration_cg04794268_eu.png"), width = 650, height = 350)
menoforest(cg04794268_eu)
dev.off() 

cg13709442 <- filter(duration_data_all, MARKER == "cg13709442")
png(here("menopause", "figures", "forest_duration_cg13709442_uk.png"), width = 650, height = 350)
menoforest(cg13709442)
dev.off() 


### timing
timing_data_path <- here("menopause", "results", "timing")

timing_results <- dir_ls(timing_data_path, regexp = "\\.txt$")

timing_data_list <- timing_results %>%
  map(~ read_tsv(
    .,
    col_names = TRUE,
    trim_ws = TRUE,
    name_repair = "unique"
  ))

# timing_data_names_raw <- names(timing_data_list)

timing_data_names <- timing_data_list %>%
  names() %>%
  str_remove_all(., paste0(timing_data_path)) %>%
  str_remove_all(., "/") %>%
  str_remove(., ".txt")

names(timing_data_list) <- timing_data_names

timing_data_list <- timing_data_list[c("timing_TWINSUK", "timing_NSHD09ken", "timing_NSHD09bell",
                                           "timing_GS1", "timing_GS2", "timing_EPICNORFOLK", "timing_FHS", 
                                           "timing_WHI", "timing_KORA", "timing_ROTTERDAM")]

timing_data_list <- lapply(timing_data_list, change_class_chr)

timing_data_all <- bind_rows(timing_data_list, .id = "STUDY")
timing_data_all$STUDY<-gsub("timing_", "", as.character(timing_data_all$STUDY))
timing_data_all$STUDY<-gsub("NSHD09ken", "NSHD2", as.character(timing_data_all$STUDY))
timing_data_all$STUDY<-gsub("NSHD09bell", "NSHD3", as.character(timing_data_all$STUDY))

timing_cpg <- c("cg04863427", "cg09366969", "cg26988423", "cg07474852")

timing_data_all <- menosubset(timing_data_all, timing_cpg)

cg04863427 <- filter(timing_data_all, MARKER == "cg04863427")
cg04863427_uk <- filter(cg04863427, STUDY %in% c("TWINSUK", "NSHD2", "NSHD3", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg04863427_us <- filter(cg04863427, STUDY %in% c("FHS", "WHI"))
cg04863427_eu <- filter(cg04863427, STUDY %in% c("KORA", "ROTTERDAM"))
menoforest(cg04863427_uk)


cg09366969 <- filter(timing_data_all, MARKER == "cg09366969")
cg09366969_uk <- filter(cg09366969, STUDY %in% c("TWINSUK", "NSHD2", "NSHD3", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg09366969_us <- filter(cg09366969, STUDY %in% c("FHS", "WHI"))
cg09366969_eu <- filter(cg09366969, STUDY %in% c("KORA", "ROTTERDAM"))
menoforest(cg09366969)
menoforest(cg09366969_uk)
menoforest(cg09366969_us)
menoforest(cg09366969_eu)

cg26988423 <- filter(timing_data_all, MARKER == "cg26988423")
cg26988423_uk <- filter(cg26988423, STUDY %in% c("TWINSUK", "NSHD2", "NSHD3", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg26988423_us <- filter(cg26988423, STUDY %in% c("FHS", "WHI"))
cg26988423_eu <- filter(cg26988423, STUDY %in% c("KORA", "ROTTERDAM"))
menoforest(cg26988423)
menoforest(cg26988423_uk)
menoforest(cg26988423_us)
menoforest(cg26988423_eu)

cg07474852 <- filter(timing_data_all, MARKER == "cg07474852")
cg07474852_uk <- filter(cg07474852, STUDY %in% c("TWINSUK", "NSHD2", "NSHD3", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg07474852_us <- filter(cg07474852, STUDY %in% c("FHS", "WHI"))
cg07474852_eu <- filter(cg07474852, STUDY %in% c("KORA", "ROTTERDAM"))
menoforest(cg07474852)
menoforest(cg07474852_uk)
menoforest(cg07474852_us)
menoforest(cg07474852_eu)

### status

status_data_names <- status_data_list %>%
  names() %>%
  str_remove_all(., paste0(status_data_path)) %>%
  str_remove_all(., "/") %>%
  str_remove(., ".txt")

names(status_data_list) <- status_data_names

status_data_list <- status_data_list[c("status_ALSPAC", "status_TWINSUK", "status_NSHD99", "status_NCDS",
                                       "status_GS1", "status_GS2", "status_EPICNORFOLK",
                                       "status_KORA", "status_ROTTERDAM")]

status_data_list <- lapply(status_data_list, change_class_chr)

status_data_all <- bind_rows(status_data_list, .id = "STUDY")
status_data_all$STUDY<-gsub("status_", "", as.character(status_data_all$STUDY))
status_data_all$STUDY<-gsub("NSHD99", "NSHD1", as.character(status_data_all$STUDY))


status_cpg <- c("cg18503679")

status_data_all <- menosubset(status_data_all, status_cpg)

cg18503679 <- filter(status_data_all, MARKER == "cg18503679")
cg18503679_uk <- filter(cg18503679, STUDY %in% c("TWINSUK", "NSHD1", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg18503679_eu <- filter(cg18503679, STUDY %in% c("KORA", "ROTTERDAM"))
menoforest(cg18503679)
menoforest(cg18503679_uk)
menoforest(cg18503679_eu)

### hrt
hrt_data_path <- here("menopause", "results", "hrt")

hrt_results <- dir_ls(hrt_data_path, regexp = "\\.txt$")

hrt_data_list <- hrt_results %>%
  map(~ read_tsv(
    .,
    col_names = TRUE,
    trim_ws = TRUE,
    name_repair = "unique"
  ))


hrt_data_names <- hrt_data_list %>%
  names() %>%
  str_remove_all(., paste0(hrt_data_path)) %>%
  str_remove_all(., "/") %>%
  str_remove(., ".txt")

names(hrt_data_list) <- hrt_data_names
 
hrt_data_list <- hrt_data_list[c("hrt_TWINSUK", "hrt_NSHD99",
                                           "hrt_GS1", "hrt_GS2", "hrt_EPICNORFOLK", "hrt_FHS", 
                                           "hrt_WHI", "hrt_KORA")]

hrt_data_list <- lapply(hrt_data_list, change_class_chr)

hrt_data_all <- bind_rows(hrt_data_list, .id = "STUDY")
hrt_data_all$STUDY<-gsub("hrt_", "", as.character(hrt_data_all$STUDY))
hrt_data_all$STUDY<-gsub("NSHD99", "NSHD1", as.character(hrt_data_all$STUDY))

hrt_cpg <- c("cg09367198", "cg06500161", "cg00049440")

hrt_data_all <- menosubset(hrt_data_all, hrt_cpg)

cg09367198 <- filter(hrt_data_all, MARKER == "cg09367198")
cg09367198_uk <- filter(cg09367198, STUDY %in% c("TWINSUK", "NSHD1", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg09367198_us <- filter(cg09367198, STUDY %in% c("FHS", "WHI"))
menoforest(cg09367198)
menoforest(cg09367198_uk)
menoforest(cg09367198_us)

cg06500161 <- filter(hrt_data_all, MARKER == "cg06500161")
cg06500161_uk <- filter(cg06500161, STUDY %in% c("TWINSUK", "NSHD1", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg06500161_us <- filter(cg06500161, STUDY %in% c("FHS", "WHI"))
menoforest(cg06500161)
menoforest(cg06500161_uk)
menoforest(cg06500161_us)

cg00049440 <- filter(hrt_data_all, MARKER == "cg00049440")
cg00049440_uk <- filter(cg00049440, STUDY %in% c("TWINSUK", "NSHD1", "NCDS", 
                                                 "GS1", "GS2", "EPICNORFOLK"))
cg00049440_us <- filter(cg00049440, STUDY %in% c("FHS", "WHI"))
menoforest(cg00049440)
menoforest(cg00049440_uk)
menoforest(cg00049440_us)