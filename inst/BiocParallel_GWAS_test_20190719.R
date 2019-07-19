# Genome-Wide BiocParallel Test
# Gabriel Odom
# 2019-07-19

######  Setup  ################################################################
library(coMethDMR)
# library(DMRcate)

###  Constant Data  ###
projectPath_char <-
  "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
pfc_df <- readRDS(
  paste0(projectPath_char, "pfc_df.RDS")
)
pfcPheno_df <- readRDS(
  paste0(projectPath_char, "pfcPheno_df.RDS")
)


###  Variable Data  ###
regions_char <- c(
  "NSHORE", "NSHELF", "SSHORE", "SSHELF", "TSS1500", "TSS200", "UTR5", "EXON1",
  "GENEBODY", "UTR3", "ISLAND"
)


######  Worker Function  ######################################################
coMethDMR_worker <- function(cpgs, beta_mat, pheno_mat, regionType_char){
  # browser()

  suppressPackageStartupMessages({
    library(coMethDMR)
    library(DMRcate)
  })


  cgi_ls <- CoMethSingleRegion(
    CpGs_char = cpgs,
    betaMatrix = beta_mat,
    method = "pearson",
    arrayType = "450k",
    returnAllCpGs = FALSE
  )
  # This second element could be a list itself
  subregions_ls <- cgi_ls[[2]]

  mods_ls <- lapply(subregions_ls, function(cgi){

    coMethBetaDF <- beta_mat[which(rownames(beta_mat) %in% cgi), ]

    lmmTest(
      betaOne_df = coMethBetaDF,
      pheno_df = pheno_mat,
      contPheno_char = "stage",
      covariates_char = c("age.brain", "sex", "Mplate", "prop.neuron"),
      modelType = "randCoef",
      arrayType = "450k"
    )

  })

  cgi_ls$modelFits_df <- do.call(rbind, mods_ls)
  cgi_ls$modelFits_df$regionType <- regionType_char
  cgi_ls

}

# Test
coMethDMR_worker(
  cpgs = c("cg20214853", "cg04677227", "cg11632906", "cg07146435"),
  beta_mat = pfc_df,
  pheno_mat = pfcPheno_df,
  regionType_char = regions_char[1]
)



######  Apply Worker over One Region Type  ####################################
region <- regions_char[1]
closeByGenomicRegion_ls <- readRDS(
  system.file(
    "extdata",
    paste0(region, "3_200.rds"),
    package = 'coMethDMR',
    mustWork = TRUE
  )
)
closeByGenomicRegion_ls <- unname(closeByGenomicRegion_ls)


library(BiocParallel)
snow_cl <- SnowParam(workers = 12, type = "SOCK")


a1 <- Sys.time()
results_ls <- bplapply(
  X = closeByGenomicRegion_ls,
  FUN = coMethDMR_worker,
  BPPARAM = snow_cl,
  beta_mat = pfc_df,
  pheno_mat = pfcPheno_df,
  regionType_char = region
)
Sys.time() - a1
# 18.8633 min over 12 cores for Islands; 7.312384 min over 12 cores for the
#   north shelf (NSHELF)



######  All Region Types  #####################################################

library(BiocParallel)
snow_cl <- SnowParam(workers = 12, type = "SOCK")

resultsAll_ls <- vector(mode = "list", length = length(regions_char))
names(resultsAll_ls) <- regions_char

aAll <- Sys.time()
for (region in regions_char) {

  print(region)

  closeByGenomicRegion_ls <- readRDS(
    system.file(
      "extdata",
      paste0(region, "3_200.rds"),
      package = 'coMethDMR',
      mustWork = TRUE
    )
  )
  closeByGenomicRegion_ls <- unname(closeByGenomicRegion_ls)

  resultsAll_ls[[region]] <- bplapply(
    X = closeByGenomicRegion_ls,
    FUN = coMethDMR_worker,
    BPPARAM = snow_cl,
    beta_mat = pfc_df,
    pheno_mat = pfcPheno_df,
    regionType_char = region
  )

}
Sys.time() - aAll
# 97.44228 min over 12 cores.

saveRDS(
  resultsAll_ls,
  file = paste0(projectPath_char, "EWAS_PFC_AllRegions_20190720.RDS")
)



######  Post-Compute Wrangling  ###############################################

projectPath_char <-
  "~/Dropbox (BBSR)/GabrielOdom/coMethDMR/vignette_parallel_computing/"
resultsAll_ls <- readRDS(
  paste0(projectPath_char, "EWAS_PFC_AllRegions_20190719.RDS")
)

# We have a giant list of lists of lists:
length(resultsAll_ls)        # 11 region types
length(resultsAll_ls$NSHORE) # 5103 North Shore Regions
resultsAll_ls$NSHORE[[1]]    # 3 elements for the first region

# Multi-layer Extraction
lmmResults_ls <- lapply(resultsAll_ls, function(regionType){

  out_ls <- lapply(regionType, `[[`, "modelFits_df")
  do.call(rbind, out_ls)

})

lmmResults_df <- do.call(rbind, lmmResults_ls)

lmmResults_df$regionType <- rep(
  names(lmmResults_ls),
  times = sapply(lmmResults_ls, nrow)
)

write.csv(
  lmmResults_df,
  file = paste0(projectPath_char, "EWAS_PFC_lmmOut_20190719.csv")
)
