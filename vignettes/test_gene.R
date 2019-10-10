data ("betasChr22_df")
data ("pheno_df")

geneList <- readRDS( system.file( "extdata", "450k_GeneByName_3.RDS", package = 'coMethDMR', mustWork = TRUE ) )

geneList$ARFGAP3

gene3_200 <- CloseBySingleRegion( CpGs_char = geneList$ARFGAP3, arrayType = "450k", maxGap = 200, minCpGs = 3 )

gene3_200

coMeth_ls <- CoMethAllRegions ( dnam = betasChr22_df, betaToM = TRUE, method = "pearson", CpGs_ls = gene3_200, arrayType = "450k",
                                returnAllCpGs = FALSE, output = "CpGs" )

coMeth_ls

results <- lmmTestAllRegions( betas = betasChr22_df,
                              region_ls = coMeth_ls,
                              pheno_df,
                              contPheno_char = "stage",
                              covariates_char = "age.brain",
                              modelType = "randCoef",
                              arrayType = "450k" )

AnnotateResults( lmmRes_df = results, arrayType = "450k" )
