
#' Fit mixed model for all regions
#'
#' @param betaMatrixAllRegions matrix of beta values for all contiguous
#'    comethylated regions, with row names = CpG ids, column names = sample ids
#' @param inFile name of input file with close by regions
#' @param fileType extension of input file: "gmt" or "RDS"
#' @param arrayType Type of array, 450k or EPIC
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param modelType model used to fit mixed model
#' @param outFile output csv file with the results for the mixed model analysis
#'
#' @return csv file with RegionID, p-value and median correlation value for
#'    each contiguous comethylated region tested.
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'    inFile = system.file("extdata", "CpGislandsChr22_ex.RDS",
#'       package = 'coMethDMR', mustWork = TRUE)
#'    lmmTestAllRegions(betaMatrixAllRegions = betaMatrixChr22_df, inFile,
#'       fileType = "RDS", arrayType = "450k", returnAllCpGs = FALSE, pheno_df,
#'       contPheno_char = "stage", covariates_char = c("age.brain", "sex"),
#'       modelType = "randCoeffMixed", outFile = "outEx.txt")
lmmTestAllRegions <- function(betaMatrixAllRegions, inFile,
                              fileType = c("gmt","RDS"),
                              arrayType = c("450k","EPIC"),
                              returnAllCpGs = FALSE, pheno_df,
                              contPheno_char, covariates_char,
                              modelType = c("randCoeffMixed", "mixed"),
                              outFile){

  ### Extract cotiguous comethylated regions ###
  coMeth_ls <- CoMethAllRegions(
    betaMatrixAllRegions, file=inFile, fileType, arrayType, returnAllCpGs
  )

  ### Create list of contiguous comethylated beta matrices ###
  CpGnames <- rownames(betaMatrixAllRegions)
  coMethBetaMatrix_ls <- lapply(
    coMeth_ls$CpGsSubregions,
    function(x) betaMatrixAllRegions[which(CpGnames%in%x), ]
  )

  ### Run mixed model for all the contiguous comethylated regions ###
  results_ls <- lapply(
    coMethBetaMatrix_ls, FUN = lmmTest, pheno_df, contPheno_char,
    covariates_char, modelType = c("randCoeffMixed", "mixed")
  )

  ### Output results ###
  out <- data.frame(
    RegionID = names(results_ls),
    Pval_Mixed_Model = unlist(unname(
      lapply(results_ls, function(x) x$`Pval_Mixed_Model`))),
    Median_Corr = unlist(unname(
      lapply(results_ls, function(x) x$`Median_Corr`)))
  )

  write.table(out, outFile, quote = FALSE)

}


