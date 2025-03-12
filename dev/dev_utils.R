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

  list(study_list, study_modified)
}

filter_studies_few_questions <- function(study_list, min_questions=10) {
  study_list |>
    purrr::keep(\(study) length(unique(study$question_id)) >= min_questions)
}

filter_study_remove_ids <- function(study_list, study_ids) {
  study_list |>
    purrr::discard(\(study) any(study$study_id %in% study_ids))
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
