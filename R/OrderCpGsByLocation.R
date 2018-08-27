
#' Title
#'
#' @param CpGs_char
#' @param CpGlocations_df
#' @param Output
#'
#' @return
#' @export
#'
#' @examples
#'
OrderCpGsByLocation <- function(CpGs_char, CpGlocations_df, Output = c("CpGs", "CpGsLocation_df")){

  ### Subset the location Data Frame  ###
  CpGs_df <- CpGlocations_df[which(CpGlocations_df$cpg%in%CpGs_char),]

  ###  Re-order this Subset  ###
  order_idx <- order(CpGs_df$CHR, as.numeric(as.character(CpGs_df$MAPINFO)))
  CpGsOrdered_df <- CpGs_df[order_idx, ]

  ### Select and return output ###
  if (Output == "CpGsLocation_df") {
    CpGsOrdered_df
    } else {
      as.character(CpGsOrdered_df$cpg)
      }


}
