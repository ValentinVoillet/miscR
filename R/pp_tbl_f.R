#' Pivot Posterior Probabilities Table
#'
#' @description
#' This function extracts the posterior mean gamma matrix from a fit object,
#' standardizes the column names by removing specific fluorochrome/channel suffixes,
#' and transforms the data into a long-format tibble.
#'
#' @details
#' The function expects `x` to be a list-like object containing a `fit` element
#' with a `mean_gamma` matrix. It uses `convert_cyt_combn_format` to standardize
#' naming before applying a regex cleanup of common flow cytometry tags.
#'
#' @param x A list object containing the model fit results (specifically `x$fit$mean_gamma`).
#' @param i Integer. The index used to extract the group name from the global or
#' environment-specific `c_obj`. Default is `1`.
#'
#' @return A \code{tibble} with four columns:
#' \item{id}{The original row names (typically observation or cell IDs)}
#' \item{combn}{The cleaned-up column names (combinations)}
#' \item{prob}{The posterior probability values}
#' \item{.grp}{The group name derived from \code{c_obj}}
#'
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom stringr str_remove_all
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
pp_tbl_f <- function(x, i = 1) {

  # Extract the mean gamma matrix
  pp_mat <- x$fit$mean_gamma

  # Standardize column names
  colnames(pp_mat) <- convert_cyt_combn_format(colnames(pp_mat), to = "std")

  # Remove fluorochrome suffixes
  patterns <- paste(
    " V450| APC| BB700| BB630| PE[-]Cy7| PE[-]CF594|",
    "BUV737| BUV395| Alx488| PE[-]Dazzle594"
  )
  colnames(pp_mat) <- stringr::str_remove_all(string = colnames(pp_mat), pattern = patterns)

  # Convert to tibble and add ID
  pp_tbl <- tibble::as_tibble(pp_mat)
  pp_tbl$id <- rownames(x$fit$mean_gamma)

  # Pivot to long format and add group metadata
  pp_tbl <- pp_tbl %>%
    tidyr::pivot_longer(
      cols = -id,
      names_to = "combn",
      values_to = "prob"
    ) %>%
    dplyr::mutate(.grp = names(c_obj)[i])

  return(pp_tbl)
}

