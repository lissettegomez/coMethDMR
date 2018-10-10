#' Find contiguous comethylated regions
#'
#' @param betaCluster_mat matrix of beta values, with rownames = sample ids,
#'    column names = CpG ids, and ordered by genomic position
#' @param rDropThresh_num min correlation between a cpg in the region with the
#'    rest of the CpGs
#'
#' @importFrom psych alpha
#'
#' @return dataframe of CpGs to keep or drop
#' @export
#'
#' @examples
#'    data(betaCluster_mat_example1)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex1)
#'    data(betaCluster_mat_example2)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex2)
#'    data(betaCluster_mat_example3)
#'    MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex3)


MarkComethylatedCpGs <- function (betaCluster_mat, rDropThresh_num = 0.5) {


  ### Calculate alpha and Store CpGs ###
  clusterAlpha_ls <- alpha(betaCluster_mat, warnings = FALSE)
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
    stringsAsFactors = FALSE
  )

  CpGs_df


}
