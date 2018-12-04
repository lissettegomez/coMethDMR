
#' Fit mixed model for all regions
#'
#' @param betaMatrixAllRegions matrix of beta values for all contiguous
#'    comethylated regions, with row names = CpG ids, column names = sample ids
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param inFile name of input file with close by regions
#' @param outFile output csv file with the results for the mixed model analysis
#' @param inFileType extension of input file: "gmt" or "RDS"
#' @param arrayType Type of array, 450k or EPIC
#' @param returnAllCpGs indicates if outputting all the CpGs in the region
#'    when there is not a contiguous comethylated region or
#'    only the CpGs in the contiguous comethylated regions
#' @param modelType model used to fit mixed model
#' @param rDropThresh_num thershold for min correlation between a cpg with sum of the
#'    rest of the CpGs
#'
#' @return text file with RegionID, p-value and median correlation value for
#'    each contiguous comethylated region tested.
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats reshape
#' @importFrom utils write.csv
#'
#' @examples
#'   data(betaMatrixChr22_df)
#'    data(pheno_df)
#'    inFile <- system.file(
#'      "extdata", "CpGislandsChr22_ex.RDS",
#'       package = 'coMethDMR', mustWork = TRUE
#'    )
#'
#'    lmmTestAllRegions(
#'      betaMatrixAllRegions = betaMatrixChr22_df,
#'      pheno_df,
#'      contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex"),
#'      inFile,
#'      inFileType = "RDS",
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE,
#'      modelType = "randCoeffMixed",
#'      rDropThresh_num = 0.5
#'    )
#'
lmmTestAllRegions <- function(betaMatrixAllRegions, pheno_df,
                              contPheno_char, covariates_char,
                              inFile, outFile = NULL,
                              inFileType = c("gmt","RDS"),
                              arrayType = c("450k","EPIC"),
                              returnAllCpGs = FALSE,
                              modelType = c("randCoeffMixed", "mixed"),
                              rDropThresh_num = 0.5){

  arrayType <- match.arg(arrayType)
  inFileType  <- match.arg(inFileType)
  modelType <- match.arg(modelType)

  ### Extract cotiguous comethylated regions ###
  coMeth_ls <- CoMethAllRegions(
    betaMatrix = betaMatrixAllRegions,
    rDropThresh_num = rDropThresh_num,
    file = inFile,
    fileType = inFileType,
    arrayType = arrayType,
    returnAllCpGs = returnAllCpGs
  )

  ### Create list of contiguous comethylated beta matrices ###
  CpGnames <- rownames(betaMatrixAllRegions)
  coMethBetaMatrix_ls <- lapply(
    coMeth_ls$CpGsSubregions,
    function(x) betaMatrixAllRegions[which(CpGnames %in% x), ]
  )

  ### Run mixed model for all the contiguous comethylated regions ###
  results_ls <- lapply(
    coMethBetaMatrix_ls,
    FUN = lmmTest,
    pheno_df, contPheno_char, covariates_char, modelType, arrayType
  )

  ### Output results ###
  outDF <- NULL
  out <- lapply(results_ls, `[[`, 1)
  for (i in 1:length(out)){
    outDF <- rbind(outDF,out[[i]])
  }

  if (is.null(outFile)){

    outDF

  } else {

    message(paste0("writing results to ", outFile))
    write.csv(outDF, outFile, quote = FALSE, row.names = FALSE)

  }

}
