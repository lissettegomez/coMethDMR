
#' Title Name a CpG region based on its genomic location
#'
#' @param CpGsOrdered_df dataframe with cpg name as character (cpg),
#' chromosome number as integer (CHR) and genomic location as integer (MAPINFO), row.names = cpg
#'
#' @return genome location of the CpGs region
#' @export
#'
#' @examples
NameRegion2 <- function(CpGsOrdered_df){

  ### Return region name based on genomic location ###
    paste0(
    "Chr", CpGsOrdered_df$CHR[1], ":",
    CpGsOrdered_df$MAPINFO[1], "-", CpGsOrdered_df$MAPINFO[length(CpGs_char)]
  )


}
