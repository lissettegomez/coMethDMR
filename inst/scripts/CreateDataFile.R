# Download and Clean GEO Data for Example
# Lissette Gomez and Lizhong Liu
# 2019-04-15
# Last edited: 2019-04-19

# Description paragraph


library(GEOquery)
library(data.table)
test <- getGEO(GEO = "GSE59685")


### Get methylation data ###
assayData <- test$GSE59685_series_matrix.txt.gz@assayData
saveRDS(assayData, "GSE59685_assay.rds")

getGEOSuppFiles(
  GEO = "GSE59685",
  baseDir = "inst/extdata",
  filter_regex = "GSE59685_betas.csv.gz"
)

betas <- fread(
  system.file ("extdata/GSE59685",
               "GSE59685_betas.csv.gz",
               package = 'coMethDMR',
               mustWork = TRUE)
)


### Get pheno data ###
pheno <- test$GSE59685_series_matrix.txt.gz@phenoData@data

pheno$stage <- as.numeric(substring(pheno$characteristics_ch1.3, 13,14))
pheno$subject.id <- as.numeric(substr(pheno$characteristics_ch1, 11, 14))
pheno$Mplate <- as.factor(substr(pheno$characteristics_ch1.1, 9,19))
pheno$sex <- as.factor (pheno$characteristics_ch1.4)
pheno$sample.id <- rownames(pheno)
pheno$age.brain <- as.numeric(substr(pheno$characteristics_ch1.6, 11, 14))
pheno$age.blood <- as.numeric(substr(pheno$characteristics_ch1.5, 11, 14))
pheno_df <- pheno[pheno$source_name_ch1 == "frontal cortex" ,
                   c("stage",  "subject.id", "Mplate", "sex", "sample.id", "age.brain")]
pheno_df <- pheno_df[!is.na(pheno_df$stage),]
saveRDS(pheno_df, "pheno_df.rds")
