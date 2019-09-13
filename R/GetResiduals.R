
#' Get Residuals
#'
#' @param beta_df data frame of beta values for all genomic regions,
#'    with row names = CpG IDs, column names = sample IDs. This is often the
#'    genome-wide array data.
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'    \code{Sample} indicating sample IDs.
#' @param covariates_char character vector for names of the covariate variables
#' @param na.action see details in na.action argument in function lm
#'
#' @return output a matrix of residual values same size as \code{beta_mat}
#'
#' @export
#'
#' @importFrom stats residuals
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(pheno_df)
#'
#'    GetResiduals(
#'      beta_df = betaMatrixChr22_df[1:10, 1:10],
#'      pheno_df = pheno_df,
#'      covariates_char = c("age.brain", "sex", "Mplate")
#'    )
#'
GetResiduals <- function(beta_df,
                         pheno_df,
                         covariates_char,
                         na.action = na.exclude){

  ### Compute M values
  mvalue_df <- log2(beta_df / (1 - beta_df))

  ### Select samples in both mvalue_df and pheno_df

  ## Make sure Sample is character but not factor
  pheno_df$Sample <- as.character(pheno_df$Sample)

  ## Check if samples in mvalue_df and pheno_df are identical
  idt <- identical(colnames(mvalue_df), pheno_df$Sample)

  if (idt){

    mvalue_df <- mvalue_df
    pheno_df <- pheno_df

  } else {

    intersectSample <- intersect(colnames(mvalue_df), pheno_df$Sample)

    ### Select samples of pheno_df based on intersect samples
    pheno_df <- pheno_df[pheno_df$Sample %in% intersectSample, ]

    ### Select samples of mvalue_df in pheno_df
    mvalue_df <- mvalue_df[ , pheno_df$Sample]

  }

  ### Create the formula
  cov_char <- paste(covariates_char, collapse = " + ")
  formula_char <- paste0("mval ~ ", cov_char)

  ### Take residuals
  resid_ls <- lapply(seq_len(nrow(mvalue_df)), function(row){

    mval <- t(mvalue_df[row, ])
    colnames(mval) <- "mval"

    dat <- cbind(mval, pheno_df)
    dat$mval <- as.numeric(dat$mval)

    fitE <- lm(formula_char, data = dat, na.action = na.action)

    residuals(fitE)

  })

  resid_df <- do.call(rbind, resid_ls)

  row.names(resid_df) <- row.names(mvalue_df)

  resid_df

}
