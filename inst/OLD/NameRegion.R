

#' Name a CpG region based on its genomic location
#'
#' @param CpGs_char vector of CpGs
#' @param CpGlocations_df dataframe with cpg name as character (cpg),
#' chromosome number as integer (CHR) and genomic location as integer (MAPINFO), row.names = cpg
#'
#' @return genome location of the CpGs region
#' @export
#'
#' @examples
#' CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#' NameRegion(CpGs_char, CpGlocations_df)
#'

NameRegion <- function(CpGs_char, CpGlocations_df=CpGlocations_df){


  ###  Subset the location Data Frame  ###
  CpGs_df <- CpGlocations_df[which(CpGlocations_df$cpg%in%CpGs_char),]

  ###  Re-order this Subset  ###
  order_idx <- order(CpGs_df$CHR, as.numeric(as.character(CpGs_df$MAPINFO)))
  CpGsOrdered_df <- CpGs_df[order_idx, ]

  ###  Return  ###
  paste0(
    "Chr", CpGsOrdered_df$CHR[1], ":",
    CpGsOrdered_df$MAPINFO[1], "-", CpGsOrdered_df$MAPINFO[length(CpGs_char)]
  )


}
