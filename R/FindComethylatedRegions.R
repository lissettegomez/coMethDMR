

#' Find contiguous comethylated regions
#'
#' @param CpGs_df a dataframe with cpg name (CpG), keep=1/drop=0 (alpha),
#'    ind=1:number of CpGs in the region
#' @param minCpGs_int n integer, minimum nubmer of cpgs for resulting regions
#'
#' @return data frame with CpG and subregion number
#' @export
#'
#' @importFrom bumphunter getSegments
#' @importFrom utils globalVariables
#'
#' @examples data(betaCluster_mat_example4)
#'    CpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaMatrix_ex4)
#'    FindComethylatedRegions(CpGs_df)
#'
FindComethylatedRegions <- function(CpGs_df, minCpGs_int = 3){

  ### Get contiguous regions of CpGs ###
  contiguousRegion_ls <- getSegments(CpGs_df$keep, cutoff = 1)
  nSegs_int <- length(contiguousRegion_ls$upIndex)

  if (nSegs_int > 0){

    ### Select segments with number of CpGs >= minCpGs ###
    contiguous_int <- lengths(contiguousRegion_ls$upIndex)
    contiguousMinCpGs_idx <- which(contiguous_int >= minCpGs_int)
    nSegsMinCpGs_int <- length(contiguousMinCpGs_idx)

    ### Create output dataframe with CpGs and contiguous comethylated subregion number ###
    #globalVariables("ind")
    ind<-NULL

       if (nSegsMinCpGs_int > 0){

         inner_ls <- lapply(1:nSegsMinCpGs_int, function(u){

           data.frame(
             CpG = subset(
               CpGs_df,
               ind %in% contiguousRegion_ls$upIndex[[contiguousMinCpGs_idx[u]]],
               select = "CpG"
             ),
             subregion = rep(
               u, length(contiguousRegion_ls$upIndex[[contiguousMinCpGs_idx[u]]])
             )
           )

         })

         contiguousRegionsCpGs <- do.call(rbind, inner_ls)


       } else {
             contiguousRegionsCpGs <- cbind(
               as.data.frame(CpGs_df$CpG),
               rep(0,length(CpGs_df$CpG))
             )
         }


  } else {
      contiguousRegionsCpGs <- cbind(
        as.data.frame(CpGs_df$CpG),
        rep(0,length(CpGs_df$CpG))
      )
    }

  colnames(contiguousRegionsCpGs) <- c("ProbeID","Subregion")

  contiguousRegionsCpGs


}
