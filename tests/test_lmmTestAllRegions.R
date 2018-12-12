# Test lmmTestAllRegions()
# Gabriel Odom
# 2018-11-28


library(coMethDMR)

data(betaMatrixChr22_df)
data(pheno_df)

CpGisland_ls <- system.file(
  "extdata", "CpGislandsChr22_ex.RDS",
  package = 'coMethDMR', mustWork = TRUE
)

coMeth_ls <- CoMethAllRegions(
  betaMatrix = betaMatrixChr22_df,
  rDropThresh_num = 0.5,
  file = CpGisland_ls,
  fileType = "RDS",
  arrayType = "450k",
  returnAllCpGs = FALSE
)

lmmTestAllRegions(
  beta_df = betaMatrixChr22_df,
  region_ls = coMeth_ls$CpGsSubregions,
  pheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex"),
  modelType = "randCoef",
  arrayType = "450k"
)
