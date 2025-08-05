library(readr)
library(dplyr)
library(moments)
library(qqman)

#####
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Locations <- as.data.frame(Locations)
Locations <- filter(Locations, chr == "chrX" | chr == "chrY")
Locations$cpg <- rownames(Locations)
write_csv(Locations[4], "../../../../Frequent_data/blacklists/450k_sex_chr.csv")
#####
sex_450 <- read_csv("../../../Frequent_data/blacklists/450k_sex_chr.csv")
sex_E <- read_csv("../../../Frequent_data/blacklists/EPIC_sex_chr_cg.csv")

input <- read_csv("../results/03_272id_NSHD99_all_resu_CurrentHRT99.csv")
#input <- input[-1]
#all(duplicated(input[,"Markername"])) #no duplicate
#g <- grep("rs",input[,"Markername"]) #no rs ID methylation probes
#w <- which(is.na(input))
input <- na.omit(input)
input <- anti_join(input, sex_E[1], by = c("Markername"="cg"))

colnames(input) <- c("MARKER","CHR","POS","BETA","SE","TSTATISTIC","PVAL","FDR","N")
write_tsv(input, "../Useful_results/A5_HRT/no_exclusioHRT_NSHD99_5A_QCed.txt")

