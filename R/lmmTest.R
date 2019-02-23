#' Fit mixed model to methylation values in one genomic region
#'
#' @param betaOne_df matrix of beta values for one genomic region,
#'    with row names = CpG IDs, column names = sample IDs
#' @param pheno_df a data frame with phenotype and covariates, with variable
#'    \code{Sample} indicating sample IDs.
#' @param contPheno_char character string of the main effect (a continuous
#'    phenotype) to be tested for association with methylation values in the
#'    region
#' @param covariates_char character vector for names of the covariate variables
#' @param modelType type of mixed model, can be \code{randCoef} for random
#'    coefficient mixed model, or \code{simple} for simple linear mixed model.
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @return  A dataframe with one row for association result of one region: \code{Estimate}, \code{StdErr}, and
#'    \code{pvalue} for the association of methylation values in the genomic
#'    region tested vs. continuous phenotype \code{contPheno_char}
#'
#' @details This function implements a mixed model to test association between
#'    methylation values in a genomic region with a continuous phenotype.
#'
#'    When \code{randCoef} is selected, the model is
#'
#'    \code{methylation M value ~ contPheno_char + covariates_char + (1|Sample) + (contPheno_char|CpG)}.
#'    The last two terms are random intercepts and slopes for each CpG.
#'
#'    When \code{simple} is selected, the model is
#'
#'    \code{methylation M value ~ contPheno_char + covariates_char + (1|Sample)}
#'
#'    In our simulation studies, we found both models are conservative, so p-values are estimated from
#'    normal distributions instead of t-distributions.
#'
#' @export
#'
#' @importFrom lmerTest lmer
#'
#' @examples
#'   data(betaMatrixChr22_df)
#'
#'   CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'
#'   coMethCpGs <- CoMethSingleRegion(CpGsChr22_char, betaMatrixChr22_df)
#'
#'   # test only the first co-methylated region
#'   coMethBetaMatrix <- betaMatrixChr22_df[coMethCpGs$CpGsSubregions[[1]], ]
#'
#'   data(pheno_df)
#'
#'   res <- lmmTest (betaOne_df = coMethBetaMatrix,
#'            pheno_df,
#'            contPheno_char = "stage",
#'            covariates_char = c("age.brain", "sex"),
#'            modelType = "randCoef",
#'            arrayType = "450k")
#'

lmmTest <- function(betaOne_df, pheno_df, contPheno_char, covariates_char,
                    modelType = c("randCoef", "simple"),
                    arrayType = c("450k","EPIC"))  {

  modelType <- match.arg(modelType)
  arrayType <- match.arg(arrayType)

  ### Transpose betaOne_df from wide to long ###
  betaOne_df$ProbeID <- row.names(betaOne_df)
  betaOneTransp_df <- reshape(
    betaOne_df,
    varying = colnames(betaOne_df)[-ncol(betaOne_df)],
    v.names = "beta",
    direction = "long",
    times = colnames(betaOne_df)[-ncol(betaOne_df)],
    timevar = "Sample"
  )

  ### Calculate M values ###
  betaOneTransp_df$Mvalue <- log2(
    betaOneTransp_df$beta / (1 - betaOneTransp_df$beta)
  )

  ### Merge transposed beta matrix with phenotype ###
  betaOnePheno_df <- merge(betaOneTransp_df, pheno_df, by = "Sample")


  ### Run the mixed model ###

  modelFormula_char <- .MakeLmmFormula(contPheno_char, covariates_char, modelType)

   tryCatch({
    f <- lmer(as.formula(modelFormula_char), betaOnePheno_df)
  }, error = function(e){ NULL })

  if(is.null(f)){

    ps_df <- data.frame(
      Estimate = NA_real_,
      StdErr = NA_real_,
      pValue = 1
    )

  } else {

    ps_mat <- coef(summary(f))[contPheno_char, c(1, 2, 4), drop = FALSE]
    ps_df <- as.data.frame(ps_mat)
    colnames(ps_df) <- c("Estimate", "StdErr", "Stat")
    rownames(ps_df) <- NULL

    ps_df$pValue <- 2 * ( 1- pnorm (abs(ps_df$Stat)))

  }

  regionName <- NameRegion(
    OrderCpGsByLocation(
      betaOne_df$ProbeID, arrayType, output = "dataframe"
    )
  )

  ### split regionName into chrom, start, end
  chrom <- sub(":.*",  "",  regionName)

  range <- sub ("c.*:", "",  regionName )

  start <- sub ("-\\d*", "", range)

  end <- sub ("\\d*.-", "", range)


  ### Return results ###

  nCpGs <- nrow(betaOne_df)

  result <- cbind (
    chrom, start, end, nCpGs,
    ps_df,
    stringsAsFactors = FALSE
  )
  result
}



.MakeLmmFormula <- function(contPheno_char, covariates_char = NULL,
                            modelType = c("randCoef", "simple")){

  modelType <- match.arg(modelType)

  baseMod_char <- "Mvalue ~ (1|Sample)"

  randomCoef_char <- paste0("(",contPheno_char, "|ProbeID)")

  if (!is.null(covariates_char)){
    cov_char <- paste(covariates_char, collapse = " + ")
    }

  ######
  if(modelType == "randCoef"){

    ifelse(
      is.null(covariates_char),
      rcMod_char <- paste(
        baseMod_char, randomCoef_char, contPheno_char, sep = " + "
      ),
      rcMod_char <- paste(
        baseMod_char, randomCoef_char, contPheno_char, cov_char, sep = " + "
      )
    )


  } else {

    ifelse(
      is.null(covariates_char),
      rcMod_char <- paste(baseMod_char, contPheno_char, sep = " + "),
      rcMod_char <- paste(baseMod_char, contPheno_char, cov_char, sep = " + ")
    )

  }


}
