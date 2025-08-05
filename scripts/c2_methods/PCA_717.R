library(readr)
info717 <- read_csv('../Info_files/717id_info_450k_ASCVD_5y.csv')
library(data.table)
meth717 <- fread('../../../Frequent_data/450k_blood_all_twins_beta_enmix_probeqced_990ids_ready.csv', data.table=FALSE)
rownames(meth717) <- meth717[, 1]
meth717 <- meth717[, -1]
meth717 <- meth717[, match(info717$Barcode, colnames(meth717))]

meth717_complete <- na.omit(meth717)

qqn <- function(i) {
  meth_qqn = qqnorm(meth717_complete[i,], plot.it=FALSE)$x
}

meth_qqn_717 <- as.data.frame(do.call(rbind, lapply(1:NROW(meth717_complete), qqn)))
rownames(meth_qqn_717) <- rownames(meth717_complete)
colnames(meth_qqn_717) <- colnames(meth717_complete)

PCA_qqn_717 <- princomp(meth_qqn_717, cor = TRUE, scores = TRUE)
PCA_qqn_717_20 <- PCA_qqn_717$loadings[, 1:20]
colnames(PCA_qqn_717_20) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
info717_20PCA <- cbind(info717, PCA_qqn_717_20)
library(ggcorrplot)
ggcorrplot(cor(info717_20PCA[, c("Age", "BMI_Nearest", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "ASCVD_100")], info717_20PCA[, 33:52]))

info239_20PCA <- read_csv('239id_info_20PCA_11062020.csv')

info_overlap_20PCA <- inner_join(info239_20PCA, info717_20PCA, by = "KCLid", suffix = c("_239", "_717"))

#PCA of 239 vs variables of 717
ggcorrplot(cor(info_overlap_20PCA[, c("Age_717", "BMI_Nearest", "CD8T_717", "CD4T_717", "NK_717", "Bcell_717", "Mono_717", "Gran_717", "ASCVD_100_717")], info_overlap_20PCA[, 34:53]))
#PCA of 717 vs variables of 239
ggcorrplot(cor(info_overlap_20PCA[, c("Age_239", "BMI", "CD8T_239", "CD4T_239", "NK_239", "Bcell_239", "Mono_239", "Gran_239", "ASCVD_100_239")], info_overlap_20PCA[, 85:104]))

