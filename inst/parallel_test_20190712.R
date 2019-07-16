# Analysis of PFC data from Lunnon et al. (2014) using Linear Mixed Model
# Lissette Gomez
# INITIAL DATE UNKNOWN
# Edits by: Gabriel Odom
# 2019-07-12

# RAM profile: from a clean system restart, R uses up to 3.9Gb, and RStudio uses
#   170Mb. After intial setup, R memory usage drops out. Total RAM usage at this
#   point is 170Mb

library(coMethDMR)
# RAM: up to 3.7Gb
library(DMRcate)
# RAM: up to 3.8Gb



######  1. import dasen normalized data  ######################################
# gitHubPath_char <-
#   "https://raw.githubusercontent.com/lissettegomez/coMethDMRPaper/master/"
# pfc_df <- readRDS(
#   url(paste0(gitHubPath_char, "data/pfc_df.RDS"))
# )
# pfcPheno_df <- readRDS(
#   url(paste0(gitHubPath_char, "data/pfcPheno_df.RDS"))
# )

# 2019-07-12; G.O.
# Lily told me that the data in the coMethDMR respository is not up-to-date. She
#   shared the updated data via DropBox.
projectPath_char <-
  "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
pfc_df <- readRDS(
  paste0(projectPath_char, "pfc_df.RDS")
)
pfcPheno_df <- readRDS(
  paste0(projectPath_char, "pfcPheno_df.RDS")
)

format(object.size(pfc_df), units = "Mb") # 441Mb
# While we have added an object that is 0.4Gb in size, the overall RAM usage
#   quickly spikes at 4.1Gb but then drops and holds at 3.6Gb.



######  2. compute co-methylated regions  #####################################
# load closeby regions that are in CGIs -> compute co-methylated regions

library(parallel)
# RAM: still at 3.6Gb


clust <- makeCluster(15)
# This spawns 15 instances of "R for Windows font-end" at 40.5Mb a piece; total
#   RAM at 608Mb + 3.6Gb = 4.2Gb
clusterEvalQ(cl = clust, library(coMethDMR))
# We now have 3.6Gb *per R instance*; total RAM at 53.8Gb + 3.6Gb = 57.4Gb. We
#   haven't even exported the data (at 0.4Gb) yet...


# G.O. - this variable is missing. It is not in any of the scripts in the
#   RScripts/ directory. I'm going to assume that Lissette meant to export the
#   two data frames.
# clusterExport(cl = clust, varlist = c("result.dir" ))
clusterExport(cl = clust, varlist = c("pfc_df", "pfcPheno_df"))
# RAM: 3.7Gb per instance; total: 55.8Gb + 3.6Gb = 59.4Gb

a <- Sys.time()
pfc_cgi_rdrop0_4_ls <- CoMethAllRegions(
  betaMatrix = pfc_df,
  regionType = "ISLAND",
  arrayType = "450k",
  rDropThresh_num = 0.4,
  returnAllCpGs = FALSE
)
Sys.time() - a
# 1.262612 hours;
#   G.O. - I'm not sure if this was over 15 cores or not. I monitored the
#   system performance, and I never saw the CPU usage break 5% (1 core out of
#   20). However, something in these last 20 lines of code is using 61Gb of
#   RAM. I'm not sure what, because the clusterExport() didn't actually export
#   anything.

saveRDS(pfc_cgi_rdrop0_4_ls, "inst/results/pfc_cgi_rdrop0_4_ls.RDS")

# R RESTART TO PURGE RAM.



######  3. Test one region  ###################################################

# Packages
library(coMethDMR)
library(DMRcate)


# Load Data and Results
pfc_cgi_ls <- readRDS("inst/results/pfc_cgi_rdrop0_4_ls.RDS")

projectPath_char <-
  "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
pfc_df <- readRDS(
  paste0(projectPath_char, "pfc_df.RDS")
)
pfcPheno_df <- readRDS(
  paste0(projectPath_char, "pfcPheno_df.RDS")
)

# Check that Mplate and sex are factors
str(pfcPheno_df)


# Compute Linear Mixed Model
oneRegion_df <- pfc_df[pfc_cgi_ls[[2]], ]
identical(as.character(pfcPheno_df$Sample), colnames(pfc_df))

a <- Sys.time()
lmm_mod <- lmmTest(
  # Data Frame of Predictor Sites
  betaOne_df = oneRegion_df,
  # Data Frame of Response and Covariates
  pheno_df = pfcPheno_df,
  # Name of the Response Column
  contPheno_char = "stage",
  # Name of the Covariate Columns
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  # Which type of model to fit?
  modelType = "randCoef",
  # Which array were sites drawn from?
  arrayType = "450k"
)
Sys.time() - a
# 1.25 seconds, but for a singular fit (region 1). I'm going to try another
#   region. For region 2, 0.5784001 sec and convergence.

# The linear mixed model output is a one-row data frame:
lmm_mod



######  4. test all regions  ##################################################
a <- Sys.time()
res_cgi_df1 <- lmmTestAllRegions(
  beta_df = pfc_df,
  region_ls = pfc_cgi_ls,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "randCoef",
  arrayType = "450k"
)
Sys.time() - a
# 34.02621 min for 4892 regions


write.csv (res_cgi_df1, "res_cgi_randCoef.csv", row.names = FALSE)

res_cgi_df2 <- lmmTestAllRegions(
  beta_df = pfc_df,
  region_ls = pfc_cgi_rdrop0_4_ls,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "simple",
  arrayType = "450k"
)

write.csv (res_cgi_df2, "res_cgi_simple.csv", row.names = FALSE)


gc()
stopCluster(clust)
