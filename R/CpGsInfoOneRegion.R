#' Test associations of individual CpGs in a genomic region with a continuous phenotype
#'
#' @param regionName_char character string of location information for a genomic region, specified in
#' the format of "chrxx:xxxxxx-xxxxxx"
#' @param betas_df data frame of beta values for all contiguous
#'    comethylated regions, with row names = CpG IDs, column names = sample IDs
#' @param pheno_df a data frame with phenotype and covariate variables, with variable "Sample" for sample IDs.
#' @param contPheno_char character string of the continuous phenotype, to be tested against methylation values
#' @param covariates_char character vector of covariate variables names
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @return a data frame with CpG ID (cpg), chromosome (chr), position (pos), and results for testing
#' association of methylation in individual CpGs with continuous phenotype (slopeEstimate, slopePval)
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom stats coef
#' @importFrom stats reshape
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'
#'    CpGsInfoOneRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      betas_df = betaMatrixChr22_df,
#'      pheno_df, contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex")
#'    )
CpGsInfoOneRegion <- function(regionName_char, betas_df, pheno_df,
                              contPheno_char, covariates_char,
                              arrayType = c("450k","EPIC")){

  arrayType <- match.arg(arrayType)

  ### Extract individual CpGs in the region ###
  CpGsToTest_char <- GetCpGsInRegion(regionName_char, arrayType = "450k")

  ### Transpose betaMatrix from wide to long ###
  CpGsBeta_df <- betas_df[
    which(rownames(betas_df) %in% CpGsToTest_char),
  ]

  ### Calculate M values ###
  CpGsMvalue_df <- log2(CpGsBeta_df / (1 - CpGsBeta_df))

  ### Match samples to test in pheno and beta data frames ###
  rownames(pheno_df) <- pheno_df$Sample
  samplesToTest <- intersect(colnames(CpGsMvalue_df), rownames(pheno_df))
  phenoTest_df <- pheno_df[samplesToTest, ]
  CpGsMvalueTest_df <- CpGsMvalue_df[ ,samplesToTest]
  #identical(rownames(phenoTest_df), colnames(CpGsMvalueTest_df))

  ### Run linear model for each CpG ###
  cov <- paste(covariates_char, collapse = "+")

  lmF <- function(Mvalue) {
    lmFormula <- as.formula(paste("Mvalue ~", contPheno_char, "+", cov ))
    tmp = coef(summary(lm(lmFormula, data=phenoTest_df)))
    tmp[contPheno_char,c(1,4)]
  }

  resultAllCpGs <- t(apply(CpGsMvalueTest_df,1,lmF))
  resultAllCpGs <- round(resultAllCpGs,4)

  ### Return results ###
  colnames(resultAllCpGs) <- c("slopeEstimate", "slopePval")
  CpGsLocation <- OrderCpGsByLocation(
    CpGs_char = CpGsToTest_char, arrayType = arrayType, output = "dataframe"
  )
  outDF <- merge(
    CpGsLocation, resultAllCpGs,
    by.x = "cpg", by.y = "row.names", sort = FALSE
  )

  outDF <- cbind(regionName_char, outDF[, -4])# remove the strand column
  colnames(outDF)[1] <- "Region"
  outDF

}


