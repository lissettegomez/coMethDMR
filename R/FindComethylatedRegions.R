


#' Find contiguous comethylated regions
#'
#' @param betaCluster_mat matrix of beta values
#' @param minCpGs_int n integer, minimum nubmer of cpgs for resulting regions
#' @param CpGlocations_df a dataframe with cpg name (cpg), chromosome (CHR) and location (MAPINFO), row.names = cpg
#' @param rDropThresh_num min correlation between a cpg in the region with the rest of the CpGs
#'
#' @importFrom bumphunter clusterMaker
#'
#' @return
#' @export
#'
#' @examples
#'


# library(psych)
# library(bumphunter)

FindComethylatedRegions <- function (betaCluster_mat,
                                  minCpGs_int=3,
                                  CpGlocations_df,
                                  rDropThresh_num=0.5) {

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

  ### Get contiguous regions of CpGs ###
  contiguousRegions <- getSegments(CpGs_df$alpha, cutoff = 1)


  contiguous_idx <- data.frame(matrix(ncol = 1, nrow = 0))

  if (length(contiguousRegions$upIndex) > 0){

    for (j in 1:length(contiguousRegions$upIndex))

    {contiguous_idx <- rbind(contiguous_idx,
                             length(contiguousRegions$upIndex[[j]]))}

    contiguosMinCpGs_idx <- as.numeric(
      rownames(subset(contiguous_idx, contiguous_idx[,1] >= minCpGs_int)))

    if (length(contiguosMinCpGs_idx) > 0){

      contiguousRegionsCpGs<-data.frame(matrix(ncol=2,nrow=0))

      for (u in 1:length(contiguosMinCpGs_idx))

      {contiguousRegionsCpGs <- rbind(
        contiguousRegionsCpGs,
        cbind(
          as.data.frame(
          subset(CpGs_df,ind %in% contiguousRegions$upIndex[[contiguosMinCpGs_idx[u]]], select="cpg")),
              rep(u, length(contiguousRegions$upIndex[[contiguosMinCpGs_idx[u]]]))))}

    } else {
      contiguousRegionsCpGs<-cbind(as.data.frame(CpGs_df$cpg),
                                   rep(0,length(CpGs_df$cpg)))}

  } else {
    contiguousRegionsCpGs<-cbind(as.data.frame(CpGs_df$cpg),
                                 rep(0,length(CpGs_df$cpg)))}

  colnames(contiguousRegionsCpGs)<-c("ProbeID","subisland")




}
