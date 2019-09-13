
#' Get Residuals
#'
#' @param dnam_df data frame of beta values for all genomic regions,
#'    with row names = CpG IDs, column names = sample IDs. This is often the
#'    genome-wide array data.
#' @param betaToM indicates if converting methylation beta values to mvalues
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'    \code{Sample} indicating sample IDs.
#' @param covariates_char character vector for names of the covariate variables
#'
#' @return output a matrix of residual values same size as \code{dnam_df}
#'
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'
#'    GetResiduals(
#'      dnam_df = betaMatrixChr22_df[1:10, 1:10],
#'      betaToM = TRUE,
#'      pheno_df = pheno_df,
#'      covariates_char = c("age.brain", "sex", "Mplate")
#'    )
#'
GetResiduals <- function(dnam_df, betaToM = TRUE,
                         pheno_df,
                         covariates_char){

  if (betaToM){
    ### Compute M values
    value_df <- log2(dnam_df / (1 - dnam_df))
  } else {
    value_df <- dnam_df
  }

  ### Select samples in both value_df and pheno_df

  ## Make sure Sample is character but not factor
  pheno_df$Sample <- as.character(pheno_df$Sample)

  ## Check if samples in value_df and pheno_df are identical
  idt <- identical(colnames(value_df), pheno_df$Sample)

  if (idt){

    value_df <- value_df
    pheno_df <- pheno_df

  } else {

    intersectSample <- intersect(colnames(value_df), pheno_df$Sample)

    ### Select samples of pheno_df based on intersect samples
    pheno_df <- pheno_df[pheno_df$Sample %in% intersectSample, ]

    ### Select samples of value_df in pheno_df
    value_df <- value_df[ , pheno_df$Sample]

  }

  ### Create the formula
  cov_char <- paste(covariates_char, collapse = " + ")
  formula_char <- paste0("val ~ ", cov_char)

  ### Take residuals
  resid_ls <- lapply(seq_len(nrow(value_df)), function(row){

    val <- t(value_df[row, ])
    colnames(val) <- "val"

    dat <- cbind(val, pheno_df)
    dat$val <- as.numeric(dat$val)

    fitE <- lm(formula_char, data = dat, na.action = na.exclude)

    residuals(fitE)

  })

  resid_df <- do.call(rbind, resid_ls)

  row.names(resid_df) <- row.names(value_df)

  resid_df

}
