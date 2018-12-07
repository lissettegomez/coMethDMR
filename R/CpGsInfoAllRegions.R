
#' Test individual CpGs in more than one region
#'
#' @param AllRegionNames_char vector of character strings with location info for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param betas_df data frame of beta values for all contiguous
#'    comethylated regions, with row names = CpG ids, column names = sample ids
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return a data frame with Region, cpg, chr, pos, slopeEstimate and slopePval for
#'    each CpG in the tested regions
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'    AllRegionNames_char <- c("chr22:18267969-18268249", "chr22:18531243-18531447")
#'
#'    CpGsInfoAllRegions(
#'      AllRegionNames_char,
#'      betas_df = betaMatrixChr22_df,
#'      pheno_df, contPheno_char = "stage",
#'      covariates_char = c("age.brain", "sex")
#'    )
CpGsInfoAllRegions <- function(AllRegionNames_char, betas_df, pheno_df,
                               contPheno_char, covariates_char,
                               arrayType = c("450k","EPIC")){

  arrayType <- match.arg(arrayType)

  resultsAllRegions_ls <- lapply(
    AllRegionNames_char,
    FUN = CpGsInfoOneRegion,
    betas_df,
    pheno_df,
    contPheno_char,
    covariates_char,
    arrayType)

  resultsAllRegions_df <- do.call(rbind, resultsAllRegions_ls)
  resultsAllRegions_df




}
