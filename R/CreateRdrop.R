#' Create R-drop values
#'
#' @param data a dataframe with rownames = sample ids, column names = CpG ids.
#' @param method correlation method, can be pearson or spearman
#'
#' @importFrom stats cor
#'
#' @return A data frame with the following columns:
#'
#' \itemize{
#'   \item{\code{CpG} : }{CpG ID}
#'   \item{\code{r_drop} : }{The correlation between each CpG with the sum of
#'   the rest of the CpGs}
#' }
#'
#' @export
#'
#' @examples
#'    data(betaMatrix_ex1)
#'    CreateRdrop(data = betaMatrix_ex1, method = "pearson")
#'
CreateRdrop <- function(data, method = c("pearson","spearman")){

  method <- match.arg(method)

  out_ls <- lapply(1:ncol(data), function(column){

    ## remove CpG i and then compute row mean
    data_no_i <- data[, -column]
    data_i <- data[, column]

    data_no_i_mean <- rowMeans(data_no_i)

    ## Correlate dat_i and dat_no_i_mean
    r.drop <- cor(data_i, data_no_i_mean, method = method)

    out_df <- data.frame(
      CpG = colnames(data)[column],
      r_drop = r.drop,
      stringsAsFactors = FALSE
    )

  })

  do.call(rbind, out_ls)

}


