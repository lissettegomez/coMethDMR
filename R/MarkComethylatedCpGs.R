


#' Find contiguous comethylated regions
#'
#' @param betaCluster_mat matrix of beta values
#' @param rDropThresh_num min correlation between a cpg in the region with the rest of the CpGs
#'
#' @importFrom psych alpha
#'
#' @return dataframe of CpGs to keep or drop
#' @export
#'
#' @examples
#'


# library(psych)
# library(bumphunter)

MarkComethylatedCpGs <- function (betaCluster_mat, rDropThresh_num=0.5) {


  ### Calculate alpha ###
  clusterAlpha_ls <- alpha(betaCluster_mat, warnings = FALSE)

  ### Drop CpGs with r.drop < threshold_r_num ###
  dropCpGs_char <- row.names(subset(clusterAlpha_ls$item.stats,
                              clusterAlpha_ls$item.stats$r.drop < rDropThresh_num))  ###drop these cpgs

  CpGs_df <- as.data.frame(rownames(clusterAlpha_ls$alpha.drop))

  colnames(CpGs_df) <- "cpg"

  CpGs_df$alpha <- ifelse(
    row.names(clusterAlpha_ls$alpha.drop) %in% dropCpGs_char, 0, 1)   ##(drop=0, keep=1)

  CpGs_df$ind <- 1:ncol(betaCluster_mat)

  CpGs_df

}
