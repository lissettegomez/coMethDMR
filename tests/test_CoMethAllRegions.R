
library(coMethDMR)
system.file(
  "extdata",
  "CpGislandsChr22_ex.RDS",
  package = 'coMethDMR',
  mustWork = TRUE
)

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
