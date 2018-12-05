
<<<<<<< HEAD
#' Fit mixed model for one region
#'
#' @param betaMatrix matrix of beta values for one contiguous comethylated region,
#'    with row names = CpG ids, column names = sample ids
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param modelType model used to fit mixed model
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return estimate, standard error and pvalue for the contiguous comethylated
#'    region being tested
=======

#' Fit mixed model to methylation values in one genomic region
#'
#' @param betaMatrix matrix of beta values for one genomic region,
#'    with row names = CpG IDs, column names = sample IDs
#'
#' @param pheno_df a data frame with phenotype and covariates, with variable \code{Sample}
#' indicating sample IDs.
#'
#' @param contPheno_char character string of the main effect (a continuous phenotype)
#' to be tested for association with methylation values in the region
#'
#' @param covariates_char character vector for names of the covariate variables
#'
#' @param modelType type of mixed model, can be \code{randCoef} for random
#' coefficient mixed model, or \code{simple} for simple linear mixed model.
#'
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @return A list with two components: (1) \code{Estimate}, \code{StdErr}, and \code{pvalue} for the association of methylation
#' values in the genomic region tested vs. continuous phenotype \code{contPheno_char};
#' (2) CpG IDs that belong to the region
#'
#' @details This function implements a mixed model to test association between methylation values in a genomic region with a continuous phenotype.
#'
#' When \code{randCoef} is selected,
#' the model is \code{methylation M value ~ contPheno_char + covariates_char + (1|Sample) + (contPheno_char|CpG)}. The last two terms are random intercepts and slopes for each CpG.
#'
#' When \code{simple} is selected, the model is \code{methylation M value ~ contPheno_char + covariates_char + (1|Sample)}
#'
>>>>>>> f6608f4c569b86a307ed9d93cd5416df47d7d2ee
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
#'   lmmTest (betaMatrix = coMethBetaMatrix,
#'            pheno_df,
#'            contPheno_char = "stage",
#'            covariates_char = c("age.brain", "sex"),
#'            modelType = "simple",
#'            arrayType = "450k")
#'

lmmTest <- function(betaMatrix, pheno_df, contPheno_char, covariates_char,
                    modelType = c("randCoef", "simple"),
                    arrayType = c("450k","EPIC"))  {

  modelType <- match.arg(modelType)
  arrayType <- match.arg(arrayType)

  ### Transpose betaMatrix from wide to long ###
  betaMatrix$ProbeID <- row.names(betaMatrix)
  betaMatrixTransp_df <- reshape(
    betaMatrix,
    varying = colnames(betaMatrix[-ncol(betaMatrix)]),
    v.names = "beta",
    direction = "long",
    times = colnames(betaMatrix[-ncol(betaMatrix)]),
    timevar = "Sample"
  )

  ### Calculate M values ###
  betaMatrixTransp_df$Mvalue <- log2(
    betaMatrixTransp_df$beta / (1 - betaMatrixTransp_df$beta)
  )

  ### Merge transposed beta matrix with phenotype ###
  betaMatrixPheno_df <- merge(betaMatrixTransp_df, pheno_df, by = "Sample")


  ### Run the mixed model ###

  modelFormula_char <- .MakeLmmFormula(contPheno_char, covariates_char, modelType)

   tryCatch({
    f <- lmer(as.formula(modelFormula_char), betaMatrixPheno_df)
  }, error = function(e){ NULL })

  if(is.null(f)){

    ps_df <- data.frame(
      Estimate = NA_real_,
      StdErr = NA_real_,
      pValue = 1
    )

  } else {

    ps_mat <- coef(summary(f))[contPheno_char, c(1, 2, 5), drop = FALSE]
    ps_df <- as.data.frame(ps_mat)
    colnames(ps_df) <- c("Estimate", "StdErr", "pValue")
    rownames(ps_df) <- NULL

  }

  regionName <- NameRegion(
    OrderCpGsByLocation(
      betaMatrix$ProbeID, arrayType, output = "dataframe"
    )
  )

  ### Return results ###
  model = cbind("Region_Name" = regionName, ps_df)
  model


}



