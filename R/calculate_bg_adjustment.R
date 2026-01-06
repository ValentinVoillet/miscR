#' Background Adjustment for Cytometry Data
#'
#' @description
#' This function identifies negative control samples (background) and subtracts their
#' percentage of positive cells from the stimulated samples. It calculates
#' `PCTPOS`, `PCTPOS_NEG`, and the adjusted value `PCTPOS_ADJ`.
#'
#' @param x A data.frame or tibble containing cytometry results.
#' @param neg_id Character. The identifier used in the `ANTIGEN` column for the negative control.
#' Default is `"negctrl"`.
#' @param group_cols A character vector of column names used to uniquely identify a biological sample
#' (e.g., Subject ID, Visit, Lab ID).
#'
#' @return A data.frame (or tibble) containing the original columns plus:
#' \item{NSUB_NEG}{Number of cells in the negative control}
#' \item{CYTNUM_NEG}{Number of positive cells in the negative control}
#' \item{PCTPOS}{Percentage of positive cells in the stimulated sample}
#' \item{PCTPOS_NEG}{Percentage of positive cells in the negative control}
#' \item{PCTPOS_ADJ}{Adjusted percentage (Stimulated - Negative Control)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' clean_data <- calculate_bg_adjustment(my_data, neg_id = "negctrl")
#' }
calculate_bg_adjustment <- function(x,
                                    neg_id = "negctrl",
                                    group_cols = c("PROTOCOL", "LABID", "ASSAYID", "PTID", "VISITNO", "RUNNUM", "SUBSET")) {

  # Ensure input is a data frame
  if (!is.data.frame(x)) stop("Input 'x' must be a data.frame or tibble.")

  # 0. Check for stimulation column (ANTIGEN or STIM)
  col_options <- c("ANTIGEN", "STIM")
  stim_col <- intersect(col_options, colnames(x))

  if (length(stim_col) == 0) {
    stop("Required column 'ANTIGEN' or 'STIM' not found in the dataset.")
  } else {
    stim_col <- stim_col[1] # Use the first match found
  }

  # 1. Extract Background (Negative Controls)
  # We create a temporary 'SAMPLE_ID' to match controls to stimulated samples
  dt_bg <- x %>%
    dplyr::filter(!!dplyr::sym(stim_col) == neg_id) %>%
    dplyr::mutate(SAMPLE_ID = do.call(paste, c(dplyr::select(., dplyr::all_of(group_cols)), sep = "_"))) %>%
    dplyr::rename(NSUB_NEG = NSUB, CYTNUM_NEG = CYTNUM) %>%
    dplyr::select(SAMPLE_ID, NSUB_NEG, CYTNUM_NEG)

  # 2. Extract Stimulated Samples
  dt_stim <- x %>%
    dplyr::filter(!!dplyr::sym(stim_col) != neg_id) %>%
    dplyr::mutate(SAMPLE_ID = do.call(paste, c(dplyr::select(., dplyr::all_of(group_cols)), sep = "_"))) %>%
    dplyr::select(
      SAMPLE_ID,
      dplyr::all_of(group_cols),
      dplyr::all_of(stim_col),
      NSUB,
      CYTNUM)

  # 3. Merge and Calculate Adjusted Percentages
  final_dt <- dt_stim %>%
    dplyr::left_join(dt_bg, by = "SAMPLE_ID") %>%
    dplyr::mutate(
      PCTPOS     = (CYTNUM / NSUB) * 100,
      PCTPOS_NEG = (CYTNUM_NEG / NSUB_NEG) * 100,
      PCTPOS_ADJ = PCTPOS - PCTPOS_NEG
    ) %>%
    dplyr::select(-SAMPLE_ID)

  return(final_dt)
}

