
#' Extract probe IDs for CpGs located in a genomic region
#'
#' @param regionName_char character string with location information for one region in
#'    this format: "chrxx:xxxxxx-xxxxxx"
#' @param arrayType Type of array, 450k or EPIC
#'
#' @return vector of CpG probe IDs mapped to the genomic region
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#'    GetCpGsInRegion(
#'      regionName_char = "chr22:18267969-18268249",
#'      arrayType = "450k"
#'    )
GetCpGsInRegion <- function(regionName_char, arrayType = c("450k","EPIC")){


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
  split1 <- unlist(strsplit(regionName_char, ":"))
  split2 <- unlist(strsplit(split1[2], "-"))
  chr <- split1[1]
  start <- split2[1]
  end <- split2[2]


  ### Subset the location Data Frame  ###
  # Remove S4 class DataFrame
  CpGlocations_df <- as.data.frame(CpGlocations_df)
  CpGlocations_df$cpg <- row.names(CpGlocations_df)
  row.names(CpGlocations_df) <- NULL
  chr_df <- CpGlocations_df[which(CpGlocations_df$chr == chr), ]
  CpGs_df <- chr_df[which(chr_df$pos >= as.integer(start) &
    chr_df$pos <= as.integer(end)), ]


  OrderCpGsByLocation(CpGs_df$cpg, arrayType = arrayType, output = "vector")

}
