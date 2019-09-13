#' Extract contiguous co-methylated genomic regions from a list of
#'   pre-defined genomic regions
#'
#' @param dnam matrix (or data frame) of beta values, with row names = CpG IDs,
#'    column names = sample IDs. This is typically genome-wide methylation beta
#'    values.
#' @param betaToM indicates if converting methylation beta values to mvalues
#' @param method method for computing correlation, can be "spearman" or "pearson"
#' @param rDropThresh_num thershold for min correlation between a cpg with sum
#'    of the rest of the CpGs
#' @param minCpGs mininum number of CpGs to be considered a "region".
#'    Only regions with more than \code{minCpGs} will be returned.
#' @param arrayType Type of array, can be "450k" or "EPIC"
#'
#' @param CpGs_ls list where each item is a character vector of CpGs IDs.
#'
#' @param regionType Type of input genomic regions (e.g. "ISLAND" for CpG island)

#' @param file an RDS or gmt file with clusters of CpG locations (i.e. CpGs
#'    located closely to each other on the genome). This file can be generated
#'    by the \code{\link{WriteCloseByAllRegions}} function.
#' @param fileType file extension for input file, can be "gmt" or "RDS"

#' @param returnAllCpGs When there is not a contiguous comethylated region in
#'    the inputing pre-defined region, \code{returnAllCpGs = 1} indicates
#'    outputting all the CpGs in the input regions, while
#'    \code{returnAllCpGs = 0} indicates not returning any CpG.
#' @param output a character vector of CpGs or a dataframe of CpGs along with
#'    rDrop info
#' @param cluster computing clusters created when using parallel computing
#'
#' @param ... Dots for internal arguments. Currently unused.
#'
#' @return  When \code{output = "dataframe"} is selected, returns a list of data frames, each with \code{CpG}
#' (CpG name), \code{Chr} (chromosome number), \code{MAPINFO} (genomic
#' position), \code{r_drop} (correlation between the CpG with rest of the
#' CpGs), \code{keep} (indicator for co-methylated CpG),
#' \code{keep_contiguous} (index for contiguous comethylated subregions).
#'
#' When \code{output = "CpGs"} is selected, returns a list, each item is a list of CpGs
#' in the contiguous co-methylated subregion.
#'
#' @details There are several ways to input genomic regions for this function: (1) use \code{CpGs_ls}
#' argument (2) use \code{regionType} argument (3) use \code{file} and \code{fileType} arguments,
#' examples of these files are at https://github.com/lissettegomez/coMethDMRdata
#'
#' @export
#'
#' @importFrom BiocParallel bplapply
#'
#' @examples
#'    data(betaMatrixChr22_df)
#'
#'
#'    CpGisland_ls <- readRDS(
#'      system.file(
#'        "extdata",
#'        "CpGislandsChr22_ex.RDS",
#'        package = 'coMethDMR',
#'        mustWork = TRUE
#'      )
#'    )
#'
#'    coMeth_ls <- CoMethAllRegions (
#'      dnam = betaMatrixChr22_df,
#'      betaToM = TRUE,
#'      method = "pearson",
#'      CpGs_ls = CpGisland_ls,
#'      arrayType = "450k",
#'      returnAllCpGs = FALSE,
#'      output = "CpGs"
#'    )
#'
#'\dontrun{
#'
#'  CoMethAllRegions (
#'    dnam = betaMatrixChr22_df,
#'    regionType = "ISLAND",
#'    arrayType = "450k",
#'    returnAllCpGs = FALSE
#'  )
#'
#'}
#'
CoMethAllRegions <- function(dnam,
                             betaToM = TRUE,
                             method = c("pearson", "spearman"),
                             rDropThresh_num = 0.4,
                             minCpGs = 3,
                             arrayType = c("450k","EPIC"),
                             CpGs_ls = NULL,
                             file = NULL,
                             fileType = c("gmt","RDS"),
                             regionType = c(
                               "ISLAND", "NSHORE", "NSHELF", "SSHORE", "SSHELF",
                               "TSS1500", "TSS200", "UTR5", "EXON1", "GENEBODY",
                               "UTR3"
                             ),
                             returnAllCpGs = FALSE,
                             output = c("CpGs", "dataframe"),
                             cluster = NULL,
                             ...){
  # browser()

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

    closeByGenomicRegion_ls <- readRDS(
      system.file(
        "extdata",
        paste0(regionType, "_3_200.rds"),
        package = 'coMethDMR',
        mustWork = TRUE
      )
    )

  }


  ### Extract contiguous comethylated region(s) from each close by region ###
  if(is.null(cluster)){

    coMethCpGsAllREgions_ls <- lapply(
      unname(closeByGenomicRegion_ls),
      FUN = CoMethSingleRegion,
      dnam = dnam,
      betaToM = betaToM,
      rDropThresh_num = rDropThresh_num,
      minCpGs = minCpGs,
      method = method,
      arrayType = arrayType,
      returnAllCpGs = returnAllCpGs
    )

  } else {

    coMethCpGsAllREgions_ls <- bplapply(
      unname(closeByGenomicRegion_ls),
      FUN = CoMethSingleRegion,
      BPPARAM = cluster,
      dnam = dnam,
      betaToM = betaToM,
      rDropThresh_num = rDropThresh_num,
      minCpGs = minCpGs,
      method = method,
      arrayType = arrayType,
      returnAllCpGs = returnAllCpGs
    )

  }


  ### return output ###
  # Return list of contiguous comethylated CpGs by Regions
  if(output == "CpGs"){

    unlist(
      lapply(coMethCpGsAllREgions_ls, `[[`, 2),
      recursive = FALSE
    )

  } else {

    out_ContigRegions <- lapply(coMethCpGsAllREgions_ls, `[[`, 1)
    out_ContigRegions[sapply(out_ContigRegions, is.null)] <- NULL
    names(out_ContigRegions) <- unlist(lapply(out_ContigRegions, `[[`, 1, 1))

    out_ContigRegions

  }


}

