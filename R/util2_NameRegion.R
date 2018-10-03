
#' Name a CpG region based on its genomic location
#'
#' @param CpGsOrdered_df dataframe with cpg name as character (cpg),
#'    chromosome number as integer (CHR) and genomic location as integer (MAPINFO)
#'
#' @return genome location of the CpGs region
#' @export
#'
#' @examples load("data/CpGsOrdered_df_ex.Rdata")
#'    NameRegion(CpGsOrdered_df)
#'
NameRegion <- function(CpGsOrdered_df){

  ### Return region name based on genomic location ###
    paste0(
    "Chr", CpGsOrdered_df$chr[1], ":",
    CpGsOrdered_df$pos[1], "-", CpGsOrdered_df$pos[nrow(CpGsOrdered_df)]
  )


}
