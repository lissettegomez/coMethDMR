# Test CoMethAllRegions()
# Gabriel Odom
# 2018-11-28

library(coMethDMR)

data(betaMatrixChr22_df)
CoMethAllRegions(
  betaMatrix = betaMatrixChr22_df,
  file = system.file(
    "extdata",
    "CpGislandsChr22_ex.RDS",
    package = 'coMethDMR',
    mustWork = TRUE),
  fileType = "RDS",
  arrayType = "450k",
  returnAllCpGs = FALSE
)
