
#' Extract individual CpGs in a region
#'
#' @param regionName_char character string with location info for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return vector of CpGs in the region
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#'    ExtractCpGs(regionName_char = "chr22:18267969-18268249", arrayType = "450k")
ExtractCpGs <- function(regionName_char, arrayType = c("450k","EPIC")){


  arrayType <- match.arg(arrayType)

  switch(arrayType,
         "450k" = {
           CpGlocations_df = IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations
         },
         "EPIC" = {
           CpGlocations_df = IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations
         }
  )


  ### Split the region name in chr and positions ###
  split1 <- strsplit(regionName_char, ":")
  split2 <- strsplit(split1[[1]][2], "-")
  split3 <- c(split1[[1]][1], split2[[1]][1], split2[[1]][2])

  ### Subset the location Data Frame  ###
  # Remove S4 class DataFrame
  CpGlocations_df <- as.data.frame(CpGlocations_df)
  CpGlocations_df$cpg <- row.names(CpGlocations_df)
  row.names(CpGlocations_df) <- NULL
  CpGs_df <- CpGlocations_df[which(
    CpGlocations_df$chr == split3[1] &
    CpGlocations_df$pos >= as.integer(split3[2]) &
    CpGlocations_df$pos <= as.integer(split3[3])),]


  OrderCpGsByLocation(CpGs_df$cpg, arrayType = arrayType, output = "vector")

}
