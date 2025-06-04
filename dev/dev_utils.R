library(CopulaSEJ)
library(energy)

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
  data_list_form <- change_value_in_study_list(studies, CopulaSEJ:::k_percentiles_to_colname(c(5, 50, 95)), 0, 0.0005)$study_list
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

change_small_value_in_study_list <- function(study_list,
                                             col_name,
                                             min_value
){
  study_modified <- vector("integer")

  for (i in seq_along(study_list)) {
    mask <- study_list[[i]][col_name] < min_value
    if (any(mask)) {
      study_modified <- append(study_modified, study_list[[i]]$study_id |> head(1))
    }
    study_list[[i]][col_name][mask] <- min_value
  }
  list(study_list=study_list, study_modified=study_modified)
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

filter_questions_general <- function(study_list, filter_fun) {
  study_list |>
    purrr::map(\(study) {
      study |> split.data.frame(study$question_id) |> purrr::keep(
        \(question_df) {
          filter_fun(question_df)
        }) |> purrr::list_rbind()
    }) |>
    purrr::keep(\(study) nrow(study) > 0)
}

filter_questions_in_studies_non_negative <- function(study_list) {
  filter_fun <- \(question_df) {
    min(question_df) >= 0
  }
  filter_questions_general(study_list, filter_fun)
}

filter_studies_support_magnitude_diff <- function(study_list, assessment_cols, max_mag_diff=1000) {
  filter_fun <- \(question_df) {
    assessments <- question_df |> dplyr::select(dplyr::all_of(assessment_cols))
    support <- c(min(assessments), max(assessments))
    diff <- support[2] / support[1]
    abs(diff) <= max_mag_diff
  }
  filter_questions_general(study_list, filter_fun)
}

filter_studies_support_diff <- function(study_list, assessment_cols, max_abs_diff=10000) {
  filter_fun <- \(question_df) {
    assessments <- question_df |> dplyr::select(dplyr::all_of(assessment_cols))
    support <- c(min(assessments), max(assessments))
    diff <- support[2] - support[1]
    abs(diff) <= max_abs_diff
  }
  filter_questions_general(study_list, filter_fun)
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

run_analysis_per_study <- function(study_list, simulation_params = NULL) {
  if (is.null(simulation_params)) {
    simulation_params <- default_simulation_params()
  }
  warnings <- list()
  results <- list()
  for (i in seq_along(study_list)) {
    study_id <- unique(study_list[[i]]$study_id) |> head(1)
    print(paste("Running study:", study_id))
    study_data <- study_list[[i]]
    study_result <- study_test_performance(study_data, simulation_params)
    warnings[[i]] <- study_result$warnings

    study_result$stats["study_id"] <- study_id
    results[[i]] <- study_result$stats
  }
  list(results = dplyr::bind_rows(results),
       warnings = warnings)
}

remove_high_magintude_difference <- function(studies, max_mag_ratio=100) {
  studies_filtered <- lapply(studies, function(study) {
    study_ranges <- study |> mutate(expert_range = `95th percentile` -`5th percentile`)
    filtered_questions <- study_ranges |> group_by(question_id) |>
      mutate(range_mag_diff = max(expert_range) / min(expert_range)) |> ungroup() |> filter(range_mag_diff < max_mag_ratio)
    return(filtered_questions)
  })
  return(studies_filtered)
}

find_high_magnitude_differences <- function(studies, max_mag_ratio = 100) {
  to_remove <- lapply(studies, function(study) {
    study_id_val <- unique(study$study_id)

    study_ranges <- study %>%
      mutate(expert_range = `95th percentile` - `5th percentile`) %>%
      group_by(question_id) %>%
      mutate(range_mag_diff = max(expert_range) / min(expert_range)) %>%
      ungroup() %>%
      filter(range_mag_diff >= max_mag_ratio) %>%
      distinct(question_id) %>%
      mutate(study_id = study_id_val)

    return(study_ranges %>% select(study_id, question_id))
  })

  # Combine all the to-be-removed entries into a single tibble
  bind_rows(to_remove)
}
