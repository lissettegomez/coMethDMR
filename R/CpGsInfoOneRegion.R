#' Test individual CpGs
#'
#' @param regionName_char character string with location info for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param betaMatrixAllRegions matrix of beta values for all contiguous
#'    comethylated regions, with row names = CpG ids, column names = sample ids
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return a data frame with cpg, chr, pos, slopeEstimate and slopePval for
#'    each CpG in the region
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
#'      betaMatrixAllRegions = betaMatrixChr22_df,
#'      pheno_df, contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex")
#'    )
CpGsInfoOneRegion <- function(regionName_char, betaMatrixAllRegions, pheno_df,
                              contPheno_char, covariates_char,
                              arrayType = c("450k","EPIC")){

  arrayType <- match.arg(arrayType)

  ### Extract individual CpGs in the region ###
  CpGsToTest_char <- CpGsInRegion(regionName_char, arrayType = "450k")

  ### Transpose betaMatrix from wide to long ###
  betaMatrixAllRegions$ProbeID <- row.names(betaMatrixAllRegions)
  CpGsBetaMatrix <- betaMatrixAllRegions[
    which(betaMatrixAllRegions$ProbeID %in% CpGsToTest_char),
  ]

  namesKept <- colnames(CpGsBetaMatrix)[-ncol(CpGsBetaMatrix)]
  CpGsBetaMatrixTransp_df <- reshape(
    CpGsBetaMatrix,
    varying = namesKept,
    v.names = "beta",
    direction = "long",
    times = namesKept,
    timevar = "Sample"
  )

  ### Calculate M values ###
  CpGsBetaMatrixTransp_df$Mvalue <- log2(
    CpGsBetaMatrixTransp_df$beta / (1 - CpGsBetaMatrixTransp_df$beta)
  )

  ### Merge transposed beta matrix with phenotype ###
  CpGsBetaMatrixPheno_df <- merge(
    CpGsBetaMatrixTransp_df, pheno_df, by = "Sample"
  )

  ### Run linear model for each CpG ###
  cov <- paste(covariates_char, collapse = "+")
  lmFormula <- as.formula(paste("Mvalue ~", contPheno_char, "+", cov ))
  resultAllCpGs <- data.frame(matrix(ncol = 3,nrow = 0))

  for (i in 1:length(CpGsToTest_char)){

    f <- lm(
      lmFormula,
      data = CpGsBetaMatrixPheno_df[
        which(CpGsBetaMatrixPheno_df$ProbeID == CpGsToTest_char[i]),
      ]
    )

    result <- coef(summary(f))[contPheno_char, c(1, 4), drop = FALSE]
    resultAllCpGs[i, ] <- cbind(CpGsToTest_char[i], round(result, 4))

  }

  ### Return results ###
  colnames(resultAllCpGs) <- c("CpG", "slopeEstimate", "slopePval")
  CpGsLocation <- OrderCpGsByLocation(
    CpGs_char = CpGsToTest_char, arrayType = arrayType, output = "dataframe"
  )
  outDF <- merge(
    CpGsLocation, resultAllCpGs,
    by.x = "cpg", by.y = "CpG", sort = FALSE
  )

  outDF[, -4] # remove the strand column

}


