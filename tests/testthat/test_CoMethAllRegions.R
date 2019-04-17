# Lizhong Liu
context("CoMethAllRegions")

data(betaMatrixChr22_df)
CpGsChr22_ls <- readRDS(
  system.file ("extdata",
               "CpGislandsChr22_ex.RDS",
               package = 'coMethDMR',
               mustWork = TRUE
  )
)

test_that("CoMethAllRegions returns df with correct classes", {

  expect_s3_class(
    CoMethAllRegions(
      betaMatrix = betaMatrixChr22_df,
      CpGs_ls = CpGsChr22_ls,
      arrayType = "450k",
      returnAllCpGs = FALSE
    ),
    "list"
  )

})
