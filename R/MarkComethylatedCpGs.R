#' Mark CpGs in contiguous and co-methylated region
#'
#' @param betaCluster_mat matrix of beta values, with rownames = sample ids,
#'    column names = CpG ids. Note that the CpGs need to be ordered by their genomic positions,
#'    this can be accomplished by the \code{OrderCpGbyLocation} function.
#'
#' @param rDropThresh_num thershold for min correlation between a cpg with sum of the
#'    rest of the CpGs
#'
#' @importFrom psych alpha
#'
#' @return A data frame with the following columns:
#'
#' \itemize{
#'   \item{\code{CpG} : }{CpG ID}
#'
#'   \item{\code{keep} : }{The CpGs with \code{keep = 1} belong to the contiguous and
#'   co-methylated region}
#'
#'   \item{\code{ind} : }{Index for the CpGs}
#'
#'   \item{\code{r_drop} : }{The correlation between each CpG with the sum of the rest of the CpGs}
#' }
#'
#' @details An outlier CpG in a genomic region will typically have low correlation with the rest of
#'  the CpGs in a genomic region. On the other hand, in a cluster of co-methylated CpGs, we expect
#'  each CpG to have high correlation with the rest of the CpGs. The \code{r.drop} statistic is used
#'  to identify these co-methylated CpGs here.
#'
#' @export
#'
#' @examples
#'    data(betaMatrix_ex1)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex1)
#'
#'    data(betaMatrix_ex2)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex2)
#'
#'    data(betaMatrix_ex3)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex3)
#'
#'    data(betaMatrix_ex4)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex4, rDropThresh_num = 0.6)
#'


MarkComethylatedCpGs <- function (betaCluster_mat, rDropThresh_num = 0.4) {

  ### Calculate alpha and Store CpGs ###

  mvalues_mat <- log2(betaCluster_mat / (1 - betaCluster_mat))

  clusterAlpha_ls <- suppressWarnings(
    alpha(mvalues_mat, warnings = FALSE)
  )

  CpGs_char <- rownames(clusterAlpha_ls$alpha.drop)

  ### Drop CpGs with r.drop < threshold_r_num ###
  # drop these cpgs
  dropCpGs_char <- row.names(
    subset(clusterAlpha_ls$item.stats,
           clusterAlpha_ls$item.stats$r.drop < rDropThresh_num)
  )

  ###  Create Output Data Frame  ###
  CpGs_df <- data.frame(
    CpG = CpGs_char,
    keep = ifelse(CpGs_char %in% dropCpGs_char, 0, 1), ##(drop=0, keep=1)
    ind = 1:ncol(betaCluster_mat),
    r_drop = clusterAlpha_ls$item.stats$r.drop,
    stringsAsFactors = FALSE
  )

  CpGs_df


}
