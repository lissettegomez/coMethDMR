#' Extract contiguous co-methylated genomic regions from a list of
#'   pre-defined genomic regions
#'
#' @param betaMatrix matrix (or data frame) of beta values, with row names = CpG IDs,
#'    column names = sample IDs. This is typically genome-wide methylation beta
#'    values.
#' @param betaToM indicates if converting methylation beta values mvalues
#' @param method method for computing correlation, can be "spearman" or "pearson"
#' @param regionType Type of input genomic regions (e.g. "ISLAND" for CpG island)
#' @param arrayType Type of array, can be "450k" or "EPIC"
#' @param file an RDS or gmt file with clusters of CpG locations (i.e. CpGs
#'    located closely to each other on the genome). This file can be generated
#'    by the \code{\link{WriteCloseByAllRegions}} function.
#' @param fileType file extension for input file, can be "gmt" or "RDS"
#' @param CpGs_ls list where each item is a character vector of CpGs IDs.
#' @param rDropThresh_num thershold for min correlation between a cpg with sum
#'    of the rest of the CpGs
#' @param returnAllCpGs When there is not a contiguous comethylated region in
#'    the inputing pre-defined region, \code{returnAllCpGs = 1} indicates
#'    outputting all the CpGs in the input regions, while
#'    \code{returnAllCpGs = 0} indicates not returning any CpG.
#' @param output a character vector of CpGs or a dataframe of CpGs along with rdrop info
#' @param ... Dots for internal arguments. Currently unused.
#'
#' @return A list of two components:
#'   \itemize{
#'     \item{\code{Contiguous_Regions} : }{A data frame with \code{CpG}
#'       (CpG name), \code{Chr} (chromosome number), \code{MAPINFO} (genomic
#'       position), \code{r_drop}(correlation between the CpG with rest of the
#'       CpGs), \code{keep} (indicator for co-methylated CpG),
#'       \code{keep_contiguous} (index for contiguous comethylated subregions)
#'     }
#'     \item{\code{CpGsSubregions} : } {results from all the regions, each item
#'       is a list of CpGs in the contiguous co-methylated subregion
#'     }
#'   }
#'
#' @export
#'
#' @examples
#'
#'
#'    data(betaMatrixChr22_df)
#'
#'    CpGsChr22_ls <- readRDS(
#'                       system.file ("extdata",
#'                                    "CpGislandsChr22_ex.RDS",
#'                                    package = 'coMethDMR',
#'                                    mustWork = TRUE
#'                                   )
#'    )
#'
#'
#'    CoMethAllRegions (
#'      betaMatrix = betaMatrixChr22_df,
#'      CpGs_ls = CpGsChr22_ls,
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE
#'    )
#'
#'
#'
#'
#'\dontrun{
#'
#'CoMethAllRegions (
#'      betaMatrix = betaMatrixChr22_df,
#'      regionType = "ISLAND",
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE
#'    )
#'
#'}
#'
CoMethAllRegions <- function(betaMatrix,
                             betaToM = TRUE,
                             regionType = c(
                               "ISLAND", "NSHORE", "NSHELF", "SSHORE", "SSHELF",
                               "TSS1500", "TSS200", "UTR5", "EXON1", "GENEBODY",
                               "UTR3"
                             ),
                             arrayType = c("450k","EPIC"),
                             file = NULL,
                             fileType = c("gmt","RDS"),
                             CpGs_ls = NULL,
                             rDropThresh_num = 0.4,
                             method = c("pearson", "spearman"),
                             returnAllCpGs = FALSE,
                             output = c("CpGs", "dataframe"),
                             ...){

  regionType <- match.arg(regionType)
  method <- match.arg(method)
  arrayType <- match.arg(arrayType)
  fileType <- match.arg(fileType)
  output <- match.arg(output)

  ### Read file of close by CpGs ###
  if(!is.null(CpGs_ls)){
    closeByGenomicRegion_ls <- CreateCpGsRegions(CpGs_ls)$regions
  } else if(!is.null(file)) {
    switch(
      fileType,
      "RDS" = {
        rdsFile <- readRDS(file)
        closeByGenomicRegion_ls <- rdsFile$regions
      },
      "gmt" = {
        gmtFile <- read_gmt(file, setType = "regions")
        closeByGenomicRegion_ls <- gmtFile$regions
      }
    )
  } else {
    closeByGenomicRegion_ls <- readRDS(system.file("extdata",
                     paste0(regionType, "3_200.rds"),
                     package = 'coMethDMR',
                     mustWork = TRUE)
    )

  }


  ### Extract contiguous comethylated region(s) from each close by region ###
  coMethCpGsAllREgions_ls <- lapply(
    unname(closeByGenomicRegion_ls),
    FUN = CoMethSingleRegion,
    betaMatrix = betaMatrix,
    betaToM = betaToM,
    rDropThresh_num = rDropThresh_num,
    method = method,
    arrayType = arrayType,
    returnAllCpGs = returnAllCpGs
  )

  ### Return list of contiguous comethylated CpGs by Regions ###
  out_ContigRegions <- lapply(coMethCpGsAllREgions_ls, `[[`, 1)
  out_ContigRegions[sapply(out_ContigRegions, is.null)] <- NULL
  names(out_ContigRegions) <- unlist(lapply(out_ContigRegions, `[[`, 1,1))

  out_coMethCpGsAll <- unlist(
    lapply(coMethCpGsAllREgions_ls, `[[`, 2),
    recursive = FALSE
  )


  ### return output ###
  if(output == "CpGs"){
    out_coMethCpGsAll
  } else {
    out_ContigRegions
  }


}

