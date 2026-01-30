#' Select Significant Marker Combinations
#'
#' @description
#' This function filters a long-format posterior probability table to identify
#' marker combinations that meet specific probability thresholds. It focuses on
#' combinations containing at least one positive marker (indicated by "+").
#'
#' @details
#' The function groups data by combination and group, then applies a quantile-based
#' filter. Specifically, it ensures that the \code{prob_min} quantile of the
#' probability distribution is at least \code{quant_min}.
#'
#' @param x A \code{tibble} or \code{data.frame} (typically the output of \code{pp_tbl_f}).
#' @param prob_min Numeric. The quantile level to evaluate (e.g., \code{0.8} for the 80th percentile).
#' Default is \code{0.8}.
#' @param quant_min Numeric. The minimum probability value required at the specified quantile.
#' Default is \code{0.25}.
#'
#' @return A character vector of unique marker combinations that passed the filtering criteria.
#'
#' @export
#'
#' @importFrom dplyr filter group_by
#' @importFrom stringr str_detect
#' @importFrom magrittr extract2
#' @importFrom stats quantile
combn_sel_vec_f <- function(x, prob_min = 0.8, quant_min = 0.25) {

  combn_sel_vec <- x %>%
    # Filter for rows containing the '+' symbol (marker presence)
    dplyr::filter(stringr::str_detect(string = combn, pattern = "\\+")) %>%
    # Group by combination and experimental group
    dplyr::group_by(combn, .grp) %>%
    # Apply quantile thresholding
    dplyr::filter(stats::quantile(prob, prob_min) >= quant_min) %>%
    # Extract the 'combn' column and return unique values
    magrittr::extract2("combn") %>%
    unique()

  return(combn_sel_vec)
}

