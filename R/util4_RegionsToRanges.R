#' Convert genomic regions in a data frame to GRanges format
#'
#' @param dat_df dataframe with variable \code{regionName} as a character vector,
#' in this format: "chrxx:xxxxxx-xxxxxx"
#'
#' @return genomic regions in GRanges format
#' @export
#'
#' @importFrom GenomicRanges GRanges
#'
#' @examples
#'  df <- data.frame ( regionName = c("chr22:19709548-19709755",
#'                                    "chr2:241721922-241722113"))
#'
#'  RegionsToRanges (df)

RegionsToRanges <- function(data_df){

  chr <- sub(":.*",  "",  data_df$regionName)

  range <- sub ("c.*:", "",  data_df$regionName )

  start <- sub ("-\\d*", "", range)

  end <- sub ("\\d*.-", "", range)

  return (GRanges ( seqname = as.factor(chr),
                    ranges = IRanges(as.numeric(start), as.numeric(end))
                  )
           )

}
