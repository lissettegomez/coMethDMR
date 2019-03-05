
#' Fit mixed model to test association between a continuous phenotype and
#' methylation values in a list of genomic regions
#'
#'
#' @param beta_df data frame of beta values for all genomic regions,
#'    with row names = CpG IDs, column names = sample IDs. This is often the
#'    genome-wide array data.
#' @param region_ls a list of genomic regions, each item is a vector of CpG IDs within a genomic region. The co-methylated
#' regions can be obtained by function \code{CoMethAllRegions}.
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'    \code{Sample} indicating sample IDs.
#' @param contPheno_char character string of the main effect (a continuous
#'    phenotype) to be tested for association with methylation values in each
#'    region
#' @param covariates_char character vector for names of the covariate variables
#' @param modelType type of mixed model, can be \code{randCoef} for random
#'    coefficient mixed model, or \code{simple} for simple linear mixed model.
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param outFile output .csv file with the results for the mixed model analysis
#'
#' @return csv file with location of the genomic region (\code{chrom, start, end}), number of CpGs (\code{nCpGs}),
#' \code{Estimate}, Standard error (\code{StdErr}) of the test statistic, p-value and False Discovery Rate (FDR)
#' for association between methylation values in each genomic region with phenotype (\code{pValue}).
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats reshape
#' @importFrom stats p.adjust
#' @importFrom utils write.csv
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'
#'    data(pheno_df)
#'
#'    CpGisland_ls <- readRDS(
#'                      system.file ("extdata",
#'                                   "CpGislandsChr22_ex.RDS",
#'                                    package = 'coMethDMR',
#'                                    mustWork = TRUE
#'                                    )
#'    )
#'
#'    coMeth_ls <- CoMethAllRegions(
#'                    betaMatrix = betaMatrixChr22_df,
#'                    file = CpGisland_ls,
#'                    fileType = "RDS",
#'                    arrayType = "450k",
#'                    rDropThresh_num = 0.4,
#'                    returnAllCpGs = FALSE
#'                )
#'
#'    lmmTestAllRegions(
#'      beta_df = betaMatrixChr22_df,
#'      region_ls = coMeth_ls$CpGsSubregions,
#'      pheno_df,
#'      contPheno_char = "stage",
#'      covariates_char = "age.brain",
#'      modelType = "randCoef",
#'      arrayType = "450k"
#'    )
#'

lmmTestAllRegions <- function(beta_df, region_ls, pheno_df,
                              contPheno_char, covariates_char,
                              modelType = c("randCoef", "simple"),
                              arrayType = c("450k","EPIC"),
                              outFile = NULL){

  modelType <- match.arg(modelType)
  arrayType <- match.arg(arrayType)

  CpGnames <- rownames(beta_df)

  coMethBetaDF_ls <- lapply(
    region_ls,
    function(x) beta_df[which(CpGnames %in% x), ]
  )

  ### Run mixed model for all the contiguous comethylated regions ###

  results_ls <- lapply(
    coMethBetaDF_ls,
    FUN = lmmTest,
    pheno_df, contPheno_char, covariates_char, modelType, arrayType
  )

  ### Output results ###

  if (length(results_ls) >0 ){

    outDF <- do.call (rbind, results_ls)

    outDF$FDR <- p.adjust(outDF$pValue, method = "fdr")


    row.names(outDF) <- NULL
  }


  if (is.null(outFile)){

    outDF

  } else {

    message(paste0("writing results to ", outFile))
    write.csv(outDF, outFile, quote = FALSE, row.names = FALSE)

  }

}
