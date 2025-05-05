library(CopulaSEJ)
library(dplyr)
library(energy)
source("dev/dev_utils.R")

calculate_decoupler_values <- function(study_data, m_func, decoupler_f, k_percentiles=c(5,50,95), flattened=TRUE) {
  if (nrow(study_data) == 0) {
    stop("Must have atleast one data row")
  }
  assessments <- study_data |> CopulaSEJ:::study_df_to_summary_array(m_func, k_percentiles)
  realizations <- CopulaSEJ:::study_df_to_realizations(study_data)
  vals <- CopulaSEJ:::assessments_to_decoupler_observations(assessments, realizations, decoupler_f)

  if (flattened) {
    return(list(decoupler_values_flat=CopulaSEJ:::flatten_3d_array_to_matrix(vals), realizations=realizations))
  } else {
    return(list(decoupler_values=vals, realizations=realizations))
  }
}

metric_list <- list(
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_CDF_decoupler()),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler()),
  list(m=CopulaSEJ:::get_median_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
  list(m=CopulaSEJ:::get_median_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler()),
  list(m=CopulaSEJ:::get_mean_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
  list(m=CopulaSEJ:::get_mean_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler())
)
output_file="dev/output/dcors_df.rds"

# run if not interactive
if (!interactive()) {
  studies <- load_data_49(relative_dev_folder=FALSE)

  metric_study_combinations <- tidyr::expand_grid(metrics=metric_list, study=studies)

  study_results <- metric_study_combinations |> purrr::pmap(\(metrics, study) {
    # calculate decoupler values
    study_id <- study$study_id |> head(1)
    values <- calculate_decoupler_values(study, metrics$m$f, metrics$decoupler$f)
    name <- glue::glue("{metrics$m$name}:{metrics$decoupler$name}")

    if (any(!is.finite(values$decoupler_values_flat))) {
      warning(glue::glue("Values are not finite for study {study_id} with metric {metrics$m$name} and decoupler {metrics$decoupler$name}. Returning NA."))
      list(
        dcor=NA,
        dcorT=NA,
        p_value_T=NA,
        name=name,
        study_id=study_id,
        dcor_test=NA,
        dcor_T_test=NA
      )
    } else {
      dcorT_result <- CopulaSEJ:::fixed_dcorT_test(values$decoupler_values_flat, values$realizations)
      dcor_result = energy::dcor.test(values$decoupler_values_flat, values$realizations, R=200)
      list(dcor=dcor_result$estimates["dCor"],
           dcorT=dcorT_result$estimate,
           p_value_T=dcorT_result$p.value,
           name=name, study_id=study_id,
           dcor_test=dcor_result,
           dcor_T_test=dcorT_result)
    }
  }) |> purrr::list_transpose()

  df_results <- tibble(!!!study_results)

  saveRDS(df_results, output_file)
}






