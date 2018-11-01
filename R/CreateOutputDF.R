

#' Create Output Dataframe
#'
#' @param keepCpGs_df a data frame with CpG = CpG name,
#'    keep = indicator for co-methylated CpG,
#'    r_drop = correlation between the CpG with rest of the CpGs
#' @param keepContiguousCpGs_df a data frame with ProbeID = CpG name,
#'    Subregion = cotiguous comethylated subregion number
#' @param CpGsOrdered_df a data frame of CpG location with
#'    chr = chromosome number, pos = genomic position, cpg = CpG name
#'
#' @return a data frame with CpG = CpG name,
#'    Chr = chromosome number,
#'    MAPINFO = genomic position,
#'    r_drop = correlation between the CpG with rest of the CpGs,
#'    keep = indicator for co-methylated CpG,
#'    keep_contiguous = cotiguous comethylated subregion number
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'    CpGsOrdered_df <- OrderCpGsByLocation(
#'    CpGsChr22_char, arrayType="450k", output = "dataframe"
#'    )
#'    betaCluster_mat <- betaMatrixChr22_df[CpGsOrdered_df$cpg,]
#'    betaClusterTransp_mat <- t(betaCluster_mat)
#'    keepCpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaClusterTransp_mat)
#'    keepContiguousCpGs_df <- FindComethylatedRegions(CpGs_df = keepCpGs_df)
#'    CreateOutputDF(keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df)
CreateOutputDF <- function(keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df){

  output_df <- merge(
    keepCpGs_df, keepContiguousCpGs_df, by.x = "CpG", by.y = "ProbeID", all.x = T
  )
  output2_df <- merge(
    CpGsOrdered_df, output_df, by.x = "cpg", by.y = "CpG"
  )

  output3_df <- output2_df[order(output2_df$chr, output2_df$pos), ]
  output3_df [is.na(output3_df)] <- 0

  coMethCpGs_df <- data.frame(
    CpG = output3_df$cpg,
    Chr = output3_df$chr,
    MAPINFO = output3_df$pos,
    r_drop = output3_df$r_drop,
    keep = output3_df$keep,
    keep_contiguous = output3_df$Subregion
  )

  coMethCpGs_df

}
