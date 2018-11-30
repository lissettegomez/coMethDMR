
#' Fit mixed model for one region
#'
#' @param betaMatrix matrix of beta values for one contiguous comethylated region,
#'    with row names = CpG ids, column names = sample ids
#' @param pheno_df a data frame with phenotype and covariates
#'    (sample ID column = "Sample")
#' @param contPheno_char character string of the phenotype name
#' @param covariates_char character vector of covariate names
#' @param modelType model used to fit mixed model
#'
#' @return list of pvalue and median correlation
#'    for the contiguous comethylated region being tested
#' @export
#'
#' @importFrom lmerTest lmer
#'
#' @examples
#'   data(betaMatrixChr22_df)
#'   CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'   coMethCpGs <- CoMethSingleRegion(CpGsChr22_char, betaMatrixChr22_df)
#'   coMethBetaMatrix <- betaMatrixChr22_df[coMethCpGs$CpGsSubregions[[1]], ]
#'   data(pheno_df)
#'   lmmTest(betaMatrix = coMethBetaMatrix, pheno_df, contPheno_char = "stage",
#'       covariates_char = c("age.brain", "sex"), modelType = "randCoeffMixed",
#'       arrayType = "450k")
lmmTest <- function(betaMatrix, pheno_df, contPheno_char, covariates_char,
                    modelType = c("randCoeffMixed", "mixed"),
                    arrayType = c("450k","EPIC"))  {

  modelType <- match.arg(modelType)

  ### Transpose betaMatrix from wide to long ###
  betaMatrix$ProbeID <- row.names(betaMatrix)
  betaMatrixTransp_df <- reshape(
    betaMatrix,
    varying = colnames(betaMatrix[-ncol(betaMatrix)]),
    v.names = "beta",
    direction = "long",
    time = colnames(betaMatrix[-ncol(betaMatrix)]),
    timevar = "Sample"
  )

  ### Calculate M values ###
  betaMatrixTransp_df$Mvalue <- log2(
    betaMatrixTransp_df$beta / (1 - betaMatrixTransp_df$beta)
  )

  ### Merge transposed beta matrix with phenotype ###
  betaMatrixPheno_df <- merge(betaMatrixTransp_df, pheno_df, by="Sample")


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

  }

  regionName <- NameRegion(
    OrderCpGsByLocation(
      betaMatrix$ProbeID, arrayType, output = "dataframe"
    )
  )

  ### Return results ###
  list(
    model = cbind("Region_Name" = regionName, ps_df),
    CpGs = betaMatrix$ProbeID
  )

}



