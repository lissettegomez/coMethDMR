

#
#' Wrapper Function to find comethyalted regions
#'
#' @param CpGs_char vector of CpGs
#' @param betaMatrix matrix of beta values, with rownames = CpG ids,
#'    column names = sample ids
#'
#' @return list of two items:
#'   1. Contiguous_Regions - a data frame of CpG location (CpG, Chr, MAPINFO),
#    r_drop (correlation between the CpG with rest of the CpGs),
#    indicator for contiguous co-methylated CpGs (keep, keep.contiguous)
#    2. CpGs_subregions - lists of CpGs in each contiguous co-methylated region
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    data(CpGsChr22_char)
#'    coMethWrapper(CpGsChr22_char, betaMatrixChr22_df)
coMethWrapper <- function(CpGs_char, betaMatrix){

  ### Order CpGs by genomic location ###
  CpGsOrdered_df <- OrderCpGsByLocation(
    CpGs_char, arrayType=c("450k"), output = "dataframe")

  ### Extract beta matrix for the input CpGs ###
  betaCluster_mat <- betaMatrix[CpGsOrdered_df$cpg,]

  ### Transpose beta matrix ###
  betaClusterTransp_mat <- t(betaCluster_mat)

  ### Mark comethylated CpGs ###
  keepCpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaClusterTransp_mat)

  ### Find contiguous comethylated regions ###
  keepContiguousCpGs_df <- FindComethylatedRegions(CpGs_df = keepCpGs_df)

  ### Split CpG dataframe by Subregion ###
  keepContiguousCpGs_ls <- SplitCpGDFbyRegion(keepContiguousCpGs_df, "450k")

  ### Create Output Data Frame  ###
  output_df <- merge(
    keepCpGs_df, keepContiguousCpGs_df, by.x = "CpG", by.y = "ProbeID", all.x = T)
  output2_df <- merge(
    CpGsOrdered_df, output_df, by.x = "cpg", by.y = "CpG")
  order_idx <- order(output2_df$chr, output2_df$pos)
  output3_df <- output2_df[order_idx,]
  output3_df [is.na(output3_df)] <- 0
  coMethCpGs_df <- data.frame(
    CpG = output3_df$cpg,
    Chr = output3_df$chr,
    MAPINFO = output3_df$pos,
    r_drop = output3_df$r_drop,
    keep = output3_df$keep,
    keep_contiguous = output3_df$Subregion)

  ### Create output list of data frame and CpGs by subregion ###
  coMethCpGs_ls <- list(coMethCpGs_df, keepContiguousCpGs_ls)
  names(coMethCpGs_ls) <- c("Contiguous_Regions", "CpGs_subregions" )

  coMethCpGs_ls

}
