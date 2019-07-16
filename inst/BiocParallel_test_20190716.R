# BiocParallel Analysis of PFC Data from Lunnon et al. (2014)
# Gabriel Odom
# 2019-07-15

# This is a replication of "inst/parallel_test_reduced_20190715.R" to implement
#   parallel computing with Snow-style clusters using BiocParallel. The first
#   script uses base parallel.

library(coMethDMR)
library(DMRcate)



######  1. import dasen normalized data  ######################################
projectPath_char <-
  "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
pfc_df <- readRDS(
  paste0(projectPath_char, "pfc_df.RDS")
)
pfcPheno_df <- readRDS(
  paste0(projectPath_char, "pfcPheno_df.RDS")
)

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

library(BiocParallel)
snow_cl <- SnowParam(workers = 12, type = "SOCK")

worker_fun <- function(x, beta_mat){

  suppressPackageStartupMessages({
    library(coMethDMR)
    library(DMRcate)
  })


  CoMethSingleRegion(
    CpGs_char = x,
    betaMatrix = beta_mat,
    method = "pearson",
    arrayType = "450k",
    returnAllCpGs = FALSE
  )

}

a <- Sys.time()
coMethCpGsAllREgions_ls <- bplapply(
  X = closeByGenomicRegion_ls,
  FUN = worker_fun,
  BPPARAM = snow_cl,
  beta_mat = pfc_df
)
Sys.time() - a
# 2.627376 min for first 480, but 1.715518 min is spent initialiazing the
#   computing cluster.
# 14.83241 min for all 4892 via BiocParallel

pfc_cgi_rdrop0_4_ls <- unlist(
  lapply(coMethCpGsAllREgions_ls, `[[`, 2),
  recursive = FALSE
)

saveRDS(pfc_cgi_rdrop0_4_ls, "inst/results/pfc_cgi_rdrop0_4_ls.RDS")
rm(pfc_cgi_rdrop0_4_ls)



######  3. lmmTest All Regions with BiocParallel  #############################

pfc_cgi_ls <- readRDS("inst/results/pfc_cgi_rdrop0_4_ls.RDS")


worker2_fun <- function(cgi, beta_mat, pheno_mat){

  suppressPackageStartupMessages({
    library(coMethDMR)
    library(DMRcate)
  })

  coMethBetaDF <- beta_mat[which(rownames(beta_mat) %in% cgi), ]

  lmmTest(
    betaOne_df = coMethBetaDF,
    pheno_df = pheno_mat,
    contPheno_char = "stage",
    covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
    modelType = "randCoef",
    arrayType = "450k"
  )

}


a1 <- Sys.time()
results_ls <- bplapply(
  X = pfc_cgi_ls,
  FUN = worker2_fun,
  BPPARAM = snow_cl,
  beta_mat = pfc_df,
  pheno_mat = pfcPheno_df
)
Sys.time() - a1
# The linear mixed model computing tops out at 4.5Gb per worker for a max RAM
#   of 58.4Gb.
# 3.735125 min over 12 cores.
# 6.800498 min over 12 cores for BiocParallel

# Post-compute wrangling
res_cgi1_df <- do.call(rbind, results_ls)
res_cgi1_df$FDR <- p.adjust(res_cgi1_df$pValue, method = "fdr")
row.names(res_cgi1_df) <- NULL

write.csv(res_cgi1_df, "inst/results/res_cgi_randCoef.csv", row.names = FALSE)
