
NameRegion <- function(CpGsOrdered_df){

  ### Return region name based on genomic location ###
    paste0(
    "Chr", CpGsOrdered_df$CHR[1], ":",
    CpGsOrdered_df$MAPINFO[1], "-", CpGsOrdered_df$MAPINFO[length(CpGs_char)]
  )


}
