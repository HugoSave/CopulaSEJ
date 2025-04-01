
load_data_49 <- function(relative_dev_folder=TRUE) {
  if (relative_dev_folder) {
    studies <- readRDS("output/data49_nov24.rds")
  } else {
    studies <- readRDS("dev/output/data49_nov24.rds")
  }
  studies
}

load_data_49_filtered <- function(relative_dev_folder=TRUE) {
  studies <- load_data_49(relative_dev_folder)
  studies <- studies |> filter_studies_few_questions(min_questions = 11)
  data_list_form <- change_value_in_study_list(studies, k_percentiles_to_colname(c(5, 50, 95)), 0, 0.001)$study_list
}

combine_lists_of_functions_to_function <- function(densities_list) {
  get_cdf_vals <- \(x) {
    purrr::map(densities_list, \(marginal) marginal$cdf(x)) |> do.call(what=cbind)
  }

  get_pdf_vals <- \(x) {
    purrr::map(densities_list, \(marginal) marginal$pdf(x)) |> do.call(what=cbind)
  }
  list(cdfs=get_cdf_vals, pdfs=get_pdf_vals)
}

#' Replace a value in a list of studies
#'
#' @param study_list list of studies.
#' @param col_name column name to replace the value in.
#' @param value_to_replace value to replace.
#' @param new_value new value to replace with.
#'
#' @return a list of studies with the value replaced and a logical vector
#' indicating indicating which studies had values replaced.
#' @export
#'
change_value_in_study_list <- function(study_list,
                                       col_name,
                                       value_to_replace,
                                       new_value) {
  study_modified <- vector("integer")

  for (i in seq_along(study_list)) {
    mask <- study_list[[i]][col_name] == value_to_replace
    if (any(mask)) {
      study_modified <- append(study_modified, study_list[[i]]$study_id |> head(1))
    }
    study_list[[i]][col_name][mask] <- new_value
  }


  list(study_list=study_list, study_modified=study_modified)
}

filter_studies_few_questions <- function(study_list, min_questions=10) {
  study_list |>
    purrr::keep(\(study) length(unique(study$question_id)) >= min_questions)
}

filter_study_remove_ids <- function(study_list, study_ids) {
  study_list |>
    purrr::discard(\(study) any(study$study_id %in% study_ids))
}

filter_study_keep_ids <- function(study_list, study_ids) {
  study_list |>
    purrr::keep(\(study) any(study$study_id %in% study_ids))
}


get_study_statistics <- function(study) {
  study_statistics <- study |>
    dplyr::summarise(
      n_questions = n_distinct(question_id),
      n_experts = n_distinct(expert_id),
      n_studies = n_distinct(study_id)
    )
  study_statistics
}
