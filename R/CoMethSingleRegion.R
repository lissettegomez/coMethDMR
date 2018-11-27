#
#' Wrapper Function to find contiguous and comethyalted sub-regions within a pre-defined genomic region
#'
#' @param CpGs_char vector of CpGs in the inputting pre-defined genomic region.
#' @param betaMatrix matrix of beta values, with row names = CpG ids,
#'    column names = sample ids
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param returnAllCpGs when there is not a contiguous comethylated region
#'
#' indicates if outputting all the CpGs in the region
#'     or
#'    only the CpGs in the contiguous comethylated regions
#' @return A list with two components:
#'
#' \itemize{
#'   \item{\code{Contiguous_Regions} : }{a data frame with \code{CpG} (CpG ID),
#'    \code{Chr} (chromosome number), \code{MAPINFO} (genomic position), \code{r_drop} (correlation between the CpG with rest of the CpGs), \code{keep} (indicator for co-methylated CpG),
#'    \code{keep_contiguous} (index for contiguous comethylated subregion)}
#'
#'   \item{\code{CpGs_subregions} : }{lists of CpGs in each contiguous co-methylated
#'       subregion}
#'
#' }
#'
#' @export
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'    CpGsChr22_char<-c("cg02953382", "cg12419862", "cg24565820", "cg04234412",
#'       "cg04824771", "cg09033563", "cg10150615", "cg18538332", "cg20007245",
#'       "cg23131131", "cg25703541")
#'    CoMethSingleRegion(CpGsChr22_char, betaMatrixChr22_df)
CoMethSingleRegion <- function(CpGs_char, betaMatrix, rDropThresh_num = 0.5,
                               arrayType = c("450k","EPIC"),
                               returnAllCpGs = TRUE){

  arrayType <- match.arg(arrayType)

  ### Order CpGs by genomic location ###
  CpGsOrdered_df <- OrderCpGsByLocation(
    CpGs_char, arrayType, output = "dataframe"
  )

  ### Extract beta matrix for the input CpGs ###
  betaCluster_mat <- betaMatrix[CpGsOrdered_df$cpg,]

  ### Transpose beta matrix ###
  betaClusterTransp_mat <- t(betaCluster_mat)

  ### Mark comethylated CpGs ###
  keepCpGs_df <- MarkComethylatedCpGs(betaCluster_mat = betaClusterTransp_mat, rDropThresh_num )

  ### Find contiguous comethylated regions ###
  keepContiguousCpGs_df <- FindComethylatedRegions(
    CpGs_df = keepCpGs_df
  )

  ### Split CpG dataframe by Subregion ###
  keepContiguousCpGs_ls <- SplitCpGDFbyRegion(
    keepContiguousCpGs_df, arrayType, returnAllCpGs
  )

  ### Create Output Data Frame  ###
  coMethCpGs_df <- CreateOutputDF(
    keepCpGs_df, keepContiguousCpGs_df, CpGsOrdered_df
  )

  ### Create output list of data frame and CpGs by subregion ###
  coMethCpGs_ls <- list(
    contiguousRegions = coMethCpGs_df,
    CpGsSubregions = keepContiguousCpGs_ls
  )


  coMethCpGs_ls

}
