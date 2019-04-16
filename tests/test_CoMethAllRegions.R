# Test CoMethAllRegions()
# Gabriel Odom
# 2018-11-28

library(coMethDMR)

data(betaMatrixChr22_df)
CpGsChr22_ls <- readRDS(
  system.file(
    "extdata",
    "CpGislandsChr22_ex.RDS",
    package = 'coMethDMR',
    mustWork = TRUE
  )
)

CoMethAllRegions(
  betaMatrix = betaMatrixChr22_df,
  CpGs_ls = CpGsChr22_ls,
  arrayType = "450k",
  returnAllCpGs = FALSE
)
