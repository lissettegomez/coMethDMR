
.MakeLmmFormula <- function(contPheno_char, covariates_char = NULL,
                            modelType = c("randCoeffMixed", "mixed")){

  modelType <- match.arg(modelType)
  baseMod_char <- "Mvalue ~ (1|Sample)"
  randomCoef_char <- paste0("(",contPheno_char, "|ProbeID)")
  cov_char <- paste(covariates_char, collapse = " + ")


  ###    ###
  if(modelType == "randCoeffMixed"){

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
