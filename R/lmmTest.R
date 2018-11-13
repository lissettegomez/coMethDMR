
#' Fit Mixed Model
#'
#' @param betaMatrix matrix of beta values for one contiguous comethylated region,
#'    with row names = CpG ids, column names = sample ids
#' @param pheno a data frame with phenotype and covariates
#' @param model model used to fit mixed model
#'
#' @return list of pvalue and median correlation
#'    for the contiguous comethylated region being tested
#' @export
#'
#' @importFrom lmerTest lmer
#'
#' @examples
lmmTest <- function(betaMatrix, pheno, method = c("randomCoefficientMixedModel", "mixedModelWithInterceptOnly"))  {

  ### Transpose betaMatrix from wide to long ###
  betaMatrix$ProbeID <- row.names(betaMatrix)
  betaMatrixTransp_df <- reshape(
    betaMatrix,
    varying = colnames(betaMatrix[-ncol(betaMatrix)]),
    v.names = "beta",
    direction = "long",
    time = colnames(betaMatrix[-ncol(betaMatrix)]),
    timevar = "Sample")

  ### Calculate M values ###
  betaMatrixTransp_df$Mvalue <- log2(betaMatrixTransp_df$beta/(1-betaMatrixTransp_df$beta))

  ### Merge transposed beta matrix with phenotype ###
  betaMatrixPheno_df <- merge(betaMatrixTransp_df, pheno, by="Sample")


  ### Run the mixed model ###
  if (model == "randomCoefficientMixedModel"){

    tryCatch({
      f<-lmer(Mvalue ~stage + (stage|ProbeID) +(1|Sample) + age.brain + sex + as.factor(Mplate), betaMatrixPheno_df)
      ps<-coef(summary(f))[2,5]
    }, error=function(e){})


  } else {


    tryCatch({
      f<-lmer(Mvalue ~stage + (1|Sample) + age.brain + sex + as.factor(Mplate), betaMatrixPheno_df)
      ps<-coef(summary(f))[2,5]
    }, error=function(e){})

  }

  ### Return results ###
  medianCorr<-median(cor(t(betaMatrix[-ncol(betaMatrix)])))
  results<-list(ps,medianCorr)
  names(results)<-c("Pval_Mixed_Model","Median_Corr")
  results


}
