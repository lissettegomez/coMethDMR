
#' Fit mixed model to test association between a continuous phenotype and
#' methylation values in a list of genomic regions
#'
#' @param beta_df data frame of beta values for all genomic regions,
#'    with row names = CpG IDs, column names = sample IDs. This is often the
#'    genome-wide array data.
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'    \code{Sample} indicating sample IDs.
#' @param contPheno_char character string of the main effect (a continuous
#'    phenotype) to be tested for association with methylation values in
#'    the region
#' @param covariates_char character vector for names of the covariate variables
#' @param modelType type of mixed model, can be \code{randCoef} for random
#'    coefficient mixed model, or \code{simple} for simple linear mixed model.
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param inFile name of input file, specifying pre-defined genomic regions
#' @param outFile output .csv file with the results for the mixed model analysis
#' @param inFileType extension of input file: "gmt" or "RDS". Sample datasets
#'    with these formats can be found in \code{inst/extdata} folder
#' @param returnAllCpGs When there is not a contiguous comethylated region in
#'    the inputing pre-defined region, \code{returnAllCpGs = 1} indicates
#'    outputting all the CpGs in the input region, \code{returnAllCpGs = 0}
#'    indicates not returning any CpG.
#' @param rDropThresh_num thershold for min correlation between a cpg with sum of the
#'    rest of the CpGs
#'
#' @return text file with RegionID, p-value for each genomic region tested.
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats reshape
#' @importFrom utils write.csv
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'    inFile <- system.file(
#'      "extdata", "CpGislandsChr22_ex.RDS",
#'       package = 'coMethDMR', mustWork = TRUE
#'    )
#'
#'    lmmTestAllRegions(
#'      beta_df = betaMatrixChr22_df,
#'      pheno_df,
#'      contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex"),
#'      inFile,
#'      inFileType = "RDS",
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE,
#'      modelType = "randCoef",
#'      rDropThresh_num = 0.5
#'    )
#'
lmmTestAllRegions <- function(beta_df, pheno_df,
                              contPheno_char, covariates_char,
                              inFile, outFile = NULL,
                              inFileType = c("gmt","RDS"),
                              arrayType = c("450k","EPIC"),
                              returnAllCpGs = FALSE,
                              modelType = c("randCoef", "simple"),
                              rDropThresh_num = 0.5){

  arrayType <- match.arg(arrayType)
  inFileType  <- match.arg(inFileType)
  modelType <- match.arg(modelType)

  ### Extract cotiguous comethylated regions ###
  coMeth_ls <- CoMethAllRegions(
    betaMatrix = beta_df,
    rDropThresh_num = rDropThresh_num,
    file = inFile,
    fileType = inFileType,
    arrayType = arrayType,
    returnAllCpGs = returnAllCpGs
  )

  ### Create list of contiguous comethylated beta matrices ###

  CpGnames <- rownames(beta_df)

  coMethBetaDF_ls <- lapply(
    coMeth_ls$CpGsSubregions,
    function(x) beta_df[which(CpGnames %in% x), ]
  )

  ### Run mixed model for all the contiguous comethylated regions ###

  results_ls <- lapply(
    coMethBetaDF_ls,
    FUN = lmmTest,
    pheno_df, contPheno_char, covariates_char, modelType, arrayType
  )

  ### Output results ###
  outDF <- NULL
  for (i in 1:length(results_ls)){
    outDF <- rbind(outDF,results_ls[[i]])
  }

  if (is.null(outFile)){

    outDF

  } else {

    message(paste0("writing results to ", outFile))
    write.csv(outDF, outFile, quote = FALSE, row.names = FALSE)

  }

}
