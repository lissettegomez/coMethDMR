# # Test lmmTestAllRegions()
# # Gabriel Odom
# # 2018-11-28
#
# # UNTIL THE OPTION TO RETURN RESULTS TO THE CONSOLE IS IMPLEMENTED, THIS TEST
# #   WILL REMAIN COMMENTED OUT.
#
# library(coMethDMR)
#
# data(betaMatrixChr22_df)
# data(pheno_df)
#
# inFile <- system.file(
#   "extdata", "CpGislandsChr22_ex.RDS",
#   package = 'coMethDMR', mustWork = TRUE
# )
#
# lmmTestAllRegions(
#   betaMatrixAllRegions = betaMatrixChr22_df,
#   pheno_df,
#   contPheno_char = "stage",
#   covariates_char = c("age.brain", "sex"),
#   inFile,
#   outFile = "outEx.txt",
#   inFileType = "RDS",
#   arrayType = "450k",
#   returnAllCpGs = FALSE,
#   modelType = "randCoeffMixed"
# )
