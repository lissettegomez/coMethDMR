# Parallel Analysis of PFC Data from Lunnon et al. (2014)
# Gabriel Odom
# 2019-07-15

# This is a modification of "inst/parallel_test.R". The original code sets up a
#   computing cluster with 15 cores, but fails to allocate these cores properly.
#   Moreover, the RAM cost to load coMethDMR is prohibitive for more than 4
#   cores. This script attempts to rectify these issues.

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


closeByGenomicRegion_ls <- readRDS(
  system.file(
    "extdata",
    "ISLAND3_200.rds",
    package = 'coMethDMR',
    mustWork = TRUE
  )
)
closeByGenomicRegion_ls <- unname(closeByGenomicRegion_ls)



######  2. compute co-methylated regions  #####################################
# load closeby regions that are in CGIs -> compute co-methylated regions

library(parallel)
# RAM: still at 3.6Gb

# clust <- makeCluster(4)
# This spawns 4 instances of "R for Windows font-end" at 40.5Mb a piece; total
#   RAM at 162Mb + 3.6Gb = 3.8Gb
# clust <- makeCluster(6)
# test for 32Gb machines
# clust <- makeCluster(8)
# This spawns 8 instances of "R for Windows font-end" at 40.5Mb a piece; total
#   RAM at 324Mb + 3.6Gb = 3.9Gb
clust <- makeCluster(12)
# total RAM: 486Mb  + 3.6Gb = 4.1Gb
clusterEvalQ(cl = clust, library(coMethDMR))
# We now have 3.6Gb *per R instance*; total RAM at 28.6Gb + 3.6Gb = 32.2Gb.
# 3.6Gb per instance + 3.6Gb master = 43Gb + 3.6Gb = 46.6Gb RAM total.
clusterEvalQ(cl = clust, library(DMRcate))
# We now have total RAM at 29Gb + 3.6Gb = 32.6Gb. We haven't even exported
#   the data (at 0.4Gb) yet...
# Up to 43.5Gb + 3.6Gb = 47.1Gb
clusterExport(
  cl = clust,
  varlist = c("pfc_df", "pfcPheno_df", "closeByGenomicRegion_ls")
)
# RAM: 4Gb per instance; total: 32.4Gb + 3.6Gb = 36Gb
# Based on this total RAM cost, I think we can go up to 8 cores (about 40Gb).
# DONE: now using 8 cores.
# For 12 cores, we are at 48.6Gb on the clusters plus 3.6Gb on the master.
#   This is 52.2Gb total, which leaves plenty of RAM available for computing.



a <- Sys.time()
# pfc_cgi_rdrop0_4_ls <- CoMethAllRegions(
#   betaMatrix = pfc_df,
#   regionType = "ISLAND",
#   arrayType = "450k",
#   returnAllCpGs = FALSE
# )
coMethCpGsAllREgions_ls <- parLapply(
  cl = clust,
  X = closeByGenomicRegion_ls,
  fun = CoMethSingleRegion,
  betaMatrix = pfc_df,
  method = "pearson",
  arrayType = "450k",
  returnAllCpGs = FALSE
)
Sys.time() - a
# 1.262612 hours serially, 91 seconds for the first 400;
# 50.33744 seconds for the first 400 over 4 cores. 1 min for first 400 over 8.
# 14.46661 min in parallel over 8 cores. My computer is tapped at 46Gb for this,
#   so I think I should try 12 cores. QUESTION: can I execute this in parallel
#   without the library loaded on the master? ANSWER: No, parLapply() can't find
#   CoMethSingleRegion(); I'll need to read up on why this is the case.

pfc_cgi_rdrop0_4_ls <- unlist(
  lapply(coMethCpGsAllREgions_ls, `[[`, 2),
  recursive = FALSE
)

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

# res_cgi_df1 <- lmmTestAllRegions(
#   beta_df = pfc_df,
#   region_ls = pfc_cgi_ls,
  # pheno_df = pfcPheno_df,
  # contPheno_char = "stage",
  # covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  # modelType = "randCoef",
  # arrayType = "450k"
# )
# # 34.02621 min for lmmTestAllRegions() over 4892 regions in serial.

a0 <- Sys.time()
CpGnames <- rownames(pfc_df)
coMethBetaDF_ls <- lapply(
  pfc_cgi_ls,
  function(x) pfc_df[which(CpGnames %in% x), ]
)
Sys.time() - a0
# 2.492514 min to split the list of data frames in serial. I'd like to
#   parallelize this too, but I think it would be worse on the memory. Exporting
#   this list is ~90Mb, while exporting the pfc_df data frame is 441Mb.

a1 <- Sys.time()
results_ls <- lapply(
  coMethBetaDF_ls,
  FUN = lmmTest,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "randCoef",
  arrayType = "450k"
)
Sys.time() - a1
# I know that if the whole computation took 35 minutes, then this probably took
#   30 minutes of it. I'm not going to test this part.

a2 <- Sys.time()
outDF <- do.call (rbind, results_ls)
outDF$FDR <- p.adjust(outDF$pValue, method = "fdr")
row.names(outDF) <- NULL
Sys.time() - a2

# RESTART R


######  5. lmmTest All Regions in Parallel  ###################################

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


# Split the PFC Data by All Known Regions
a0 <- Sys.time()
CpGnames <- rownames(pfc_df)
coMethBetaDF_ls <- lapply(
  pfc_cgi_ls,
  function(x) pfc_df[which(CpGnames %in% x), ]
)
Sys.time() - a0
# 2.492514 min in serial. Don't distribute this.


# Set up cluster
library(parallel)
clust <- makeCluster(12)
clusterEvalQ(cl = clust, library(coMethDMR))
clusterEvalQ(cl = clust, library(DMRcate))
clusterExport(
  cl = clust,
  varlist = c("coMethBetaDF_ls", "pfcPheno_df")
)
# 3.6Gb over 12 clusters + 4.4Gb master = 43.5Gb + 4.4Gb = 47.9Gb


# Parallel LMM
a1 <- Sys.time()
results_ls <- parLapply(
  cl = clust,
  coMethBetaDF_ls,
  fun = lmmTest,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "randCoef",
  arrayType = "450k"
)
Sys.time() - a1
# The linear mixed model computing tops out at 4.5Gb per worker for a max RAM
#   of 58.4Gb.
# 3.735125 min over 12 cores.


# Post-compute wrangling
res_cgi1_df <- do.call(rbind, results_ls)
res_cgi1_df$FDR <- p.adjust(res_cgi1_df$pValue, method = "fdr")
row.names(res_cgi1_df) <- NULL

write.csv(res_cgi1_df, "inst/results/res_cgi_randCoef.csv", row.names = FALSE)
rm(results_ls, res_cgi1_df)

# If you also care about the simple linear mixed model (instead of the random
#   coefficient version), run that while the clusters are still active and have
#   all the data you already need on them.
a2 <- Sys.time()
results2_ls <- parLapply(
  cl = clust,
  coMethBetaDF_ls,
  fun = lmmTest,
  pheno_df = pfcPheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
  modelType = "simple",
  arrayType = "450k"
)
Sys.time() - a2
# 2.473109 min over 12 cores.


# Post-compute wrangling
res_cgi2_df <- do.call(rbind, results2_ls)
res_cgi2_df$FDR <- p.adjust(res_cgi2_df$pValue, method = "fdr")
row.names(res_cgi2_df) <- NULL

write.csv (res_cgi2_df, "inst/results/res_cgi_simple.csv", row.names = FALSE)
rm(results2_ls, res_cgi2_df)

stopCluster(clust)
# Now we are back to 4Gb of RAM in the single R+RStudio session.
