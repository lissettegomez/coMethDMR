
#' Create a list object with class CpGsRegions
#'
#' @param CpGs_ls a list where each item is a character vector of CpGs IDs in a
#'    region. Each vector should be named with the region name.
#'
#' @return A list object with class CpGsRegions
#' @export
#'
#' @examples
#'
#'    CpGsChr22_ls <- readRDS(
#'                       system.file ("extdata",
#'                                    "CpGislandsChr22_ex.RDS",
#'                                    package = 'coMethDMR',
#'                                    mustWork = TRUE
#'                                   )
#'    )
#'    CreateCpGsRegions(CpGsChr22_ls)
#'
CreateCpGsRegions <- function(CpGs_ls){
  out_PC <- CreatePathwayCollection(
    sets_ls = CpGs_ls,
    TERMS = names(CpGs_ls),
    setType = "regions"
  )
  classOut <- class(out_PC)
  class(out_PC) <- c("CpGsRegions", classOut)
  out_PC
}
