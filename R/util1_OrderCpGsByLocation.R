
#' Order CpGs by genomic location
#'
#' @param CpGs_char vector of CpGs
#' @param arrayType Type of array, 450k or EPIC
#' @param Output vector of CpGs or dataframe with CpGs, CHR, MAPINFO
#'
#' @return vector of ordered CpGs or dataframe with ordered CpGs, CHR, MAPINFO
#' @export
#'
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#'
#' @examples
#' CpGs_char <- c("cg04677227", "cg07146435", "cg11632906", "cg20214853")
#' OrderCpGsByLocation(CpGs_char, arrayType=c("EPIC"), output = "dataframe")
#'
OrderCpGsByLocation <- function(CpGs_char, arrayType=c("450k","EPIC"), output = c("vector", "dataframe")){

  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = data(IlluminaHumanMethylation450kanno.ilmn12.hg19),
         "EPIC" = data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))

  ### Subset the location Data Frame  ###
  #data(Locations)
  CpGlocations_df <- Locations
  CpGlocations_df$cpg <- row.names(CpGlocations_df)
  row.names(CpGlocations_df) <- NULL
  CpGs_df <- CpGlocations_df[which(CpGlocations_df$cpg%in%CpGs_char),]

  ###  Re-order this Subset  ###
  order_idx <- order(CpGs_df$chr, CpGs_df$pos)
  CpGsOrdered_df <- CpGs_df[order_idx, ]

  ### Select and return output ###
  if (output == "dataframe") {
    CpGsOrdered_df
  } else {
    as.character(CpGsOrdered_df$cpg)
  }


}
