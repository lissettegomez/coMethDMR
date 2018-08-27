
#' Order CpGs by genomic location
#'
#' @param CpGs_char vector of CpGs
#' @param CpGlocations_df dataframe with cpg name as character (cpg),
#' chromosome number as integer (CHR) and genomic location as integer (MAPINFO), row.names = cpg
#' @param vector Output: vector of CpGs or dataframe with CpGs, CHR, MAPINFO
#'
#' @return of ordered CpGs or dataframe of ordered CpGs, CHR, MAPINFO
#' @export
#'
#' @examples
#'
OrderCpGsByLocation <- function(CpGs_char, CpGlocations_df, Output = c("CpGs", "CpGsOrdered_df")){

  ### Subset the location Data Frame  ###
  CpGs_df <- CpGlocations_df[which(CpGlocations_df$cpg%in%CpGs_char),]

  ###  Re-order this Subset  ###
  order_idx <- order(CpGs_df$CHR, as.numeric(as.character(CpGs_df$MAPINFO)))
  CpGsOrdered_df <- CpGs_df[order_idx, ]

  ### Select and return output ###
  if (Output == "CpGsOrdered_df") {
    CpGsOrdered_df
    } else {
      as.character(CpGsOrdered_df$cpg)
      }


}
