# Test lmmTestAllRegions()
# Gabriel Odom
# 2018-11-28

data(betaMatrixChr22_df)
data(pheno_df)

inFile <- system.file(
  "extdata", "CpGislandsChr22_ex.RDS",
  package = 'coMethDMR', mustWork = TRUE
)

lmmTestAllRegions(
  betaMatrixAllRegions = betaMatrixChr22_df,
  pheno_df,
  contPheno_char = "stage",
  covariates_char = c("age.brain", "sex"),
  inFile,
  outFile = "outEx.txt",
  inFileType = "RDS",
  arrayType = "450k",
  returnAllCpGs = FALSE,
  modelType = "randCoeffMixed"
)
