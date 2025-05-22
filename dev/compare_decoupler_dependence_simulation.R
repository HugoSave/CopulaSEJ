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
    return(list(decoupler_values_flat=CopulaSEJ:::flatten_3d_array_to_matrix(vals),
                realizations=realizations,
                assessments_flat=CopulaSEJ:::flatten_3d_array_to_matrix(assessments)
                ))
  } else {
    return(list(decoupler_values=vals,
                realizations=realizations,
                assessments=assessments
                ))
  }
}

metric_list <- list(
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_CDF_decoupler()),
  #list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_linear_decoupler(D_tilde = 3, compose_sigmoid = FALSE, short_name="Linear")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_linear_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_E", short_name="Linear_MnE")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_linear_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_G", short_name="Linear_MnG")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_linear_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "median", short_name="Linear_Md")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_relative_decoupler(D_tilde = 3, compose_sigmoid = FALSE, short_name="Relative")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_relative_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_E", short_name="Relative_MnE")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_relative_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_G", short_name="Relative_MnG")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_relative_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "median", short_name="Relative_Md")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_support_ratio_decoupler(global_support = TRUE)),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_ratio_decoupler(D_tilde = 3, compose_sigmoid = FALSE, short_name="Ratio")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_ratio_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "median", short_name="Ratio_Md")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_ratio_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_E", short_name="Ratio_MnE")),
  list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_ratio_decoupler(D_tilde = 1, compose_sigmoid = FALSE, m_preprocess = "mean_G", short_name="Ratio_MnG"))
  # list(m=CopulaSEJ:::get_mean_summarizing_function(unified_support=TRUE), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
  # list(m=CopulaSEJ:::get_mean_summarizing_function(unified_support=TRUE), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler())
)

# metric_list <- list(
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_CDF_decoupler()),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler()),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_mean_linear_decoupler(global_support = TRUE)),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_mean_linear_decoupler(global_support=FALSE)),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_mean_ratio_decoupler(global_support = TRUE)),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_mean_ratio_decoupler(global_support = FALSE)),
#   list(m=CopulaSEJ:::get_three_quantiles_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_support_ratio_decoupler(global_support = TRUE)),
#   list(m=CopulaSEJ:::get_median_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
#   list(m=CopulaSEJ:::get_median_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler()),
#   list(m=CopulaSEJ:::get_median_summarizing_function(), decoupler=CopulaSEJ:::get_sigmoid_relative_decoupler())
#   # list(m=CopulaSEJ:::get_mean_summarizing_function(unified_support=TRUE), decoupler=CopulaSEJ:::get_sigmoid_linear_decoupler()),
#   # list(m=CopulaSEJ:::get_mean_summarizing_function(unified_support=TRUE), decoupler=CopulaSEJ:::get_sigmoid_ratio_decoupler())
# )

run_decoupler_summarizer_dependence_analysis <- function(
  studies,
  metric_list,
  output_file,
  reject_experts=FALSE,
  rejection_threshold=0.05,
  rejection_test="distance_correlation",
  min_nr_experts=0
) {
  # check if the output file already exists

  metric_study_combinations <- tidyr::expand_grid(metrics=metric_list, study=studies)

  study_results <- metric_study_combinations |> purrr::pmap(\(metrics, study) {
    # calculate decoupler values
    study_id <- study$study_id |> head(1)
    nr_experts <- study$expert_id |> unique() |> length()
    nr_questions <- study$question_id |> unique() |> length()
    realizations <- matrix(study_df_to_realizations(study), nrow=nr_questions, ncol=1)
    experts_rejected = NA
    nr_experts_rejected = NA
    accepted_experts = seq_len(nr_experts)
    values <- calculate_decoupler_values(study, metrics$m$f, metrics$decoupler$f, flattened=FALSE)
    if (reject_experts) {
      rejection_results <- reject_experts(values$assessments, values$realizations, rejection_threshold, test=rejection_test,
                                          decoupler=metrics$decoupler, min_nr_experts=min_nr_experts)
      if (length(rejection_results$accepted_experts) == 0) {
        stop("No experts accepted. Please check the rejection threshold or set rejection_min_experts > 0.")
      }
      accepted_experts = rejection_results$accepted_experts
      experts_rejected = rejection_results$rejected_experts
      nr_experts_rejected = length(rejection_results$rejected_experts)
    }
    acceped_Z <- values$decoupler_values[,accepted_experts,,drop=FALSE]
    Z_flat <- CopulaSEJ:::flatten_3d_array_to_matrix(acceped_Z)
    accepted_M <- values$assessments[,accepted_experts,,drop=FALSE]
    M_flat <- CopulaSEJ:::flatten_3d_array_to_matrix(accepted_M)
    name <- glue::glue("{metrics$m$name}:{metrics$decoupler$name}")
    D=dim(values$assessments)[3]
    D_tilde=dim(values$decoupler_values)[3]

    if (any(!is.finite(Z_flat))) {
      warning(glue::glue("Values are not finite for study {study_id} with metric {metrics$m$name} and decoupler {metrics$decoupler$name}. Returning NA."))
      tibble::tibble_row(
        dcor=NA,
        dcorT=NA,
        p_value_T=NA,
        E=nr_experts,
        D=D,
        D_tilde=D_tilde,
        name=name,
        m=metrics$m$short_name,
        decoupler=metrics$decoupler$short_name,
        study_id=study_id,
        dcor_test=NA,
        dcor_T_test=NA,
        experts_rejected=list(experts_rejected),
        nr_experts_rejected=nr_experts_rejected,
        nr_questions=nr_questions
      )
    } else {
      dcorT_test_M <- CopulaSEJ:::fixed_dcorT_test(Z_flat, M_flat)
      dcorT_test_Q <- CopulaSEJ:::fixed_dcorT_test(Z_flat, realizations)
      dcor_result = energy::dcor.test(Z_flat, M_flat, R=200)
      tibble::tibble_row(dcor_M=dcor_result$estimates["dCor"],
           dcorT_M=dcorT_test_M$estimate,
           dcorT_Q=dcorT_test_Q$estimate,
           p_value_M=dcorT_test_M$p.value,
           p_value_Q=dcorT_test_Q$p.value,
           E=nr_experts,
           D=D,
           D_tilde=D_tilde,
           name=name, study_id=study_id,
           m=metrics$m$short_name,
           decoupler=metrics$decoupler$short_name,
           dcor_test_M=dcor_result,
           dcor_T_test_M=dcorT_test_M,
           dcorT_test_Q=dcorT_test_Q,
           experts_rejected=list(experts_rejected),
           nr_experts_rejected=nr_experts_rejected,
           nr_questions=nr_questions
      )
    }
  }, .progress = TRUE) |> purrr::list_rbind()

  #browsr
  #df_results <- tibble(!!!study_results)

  print(glue::glue("Saving results to {output_file}"))
  saveRDS(study_results, output_file)
}

decoupler_test_output_file_name <- function(
  reject_experts=FALSE,
  rejection_threshold=0.05,
  date=Sys.Date(),
  rel_dev=FALSE
) {
  rej_expert_str <- ifelse(reject_experts, paste0("RejE(", rejection_threshold, ")"), "noR")
  dev_prefix <- ifelse(rel_dev, "", "dev/")
  glue::glue("{dev_prefix}output/dependence_test_{rej_expert_str}_{date}.rds")
}

# run if not interactive
if (!interactive()) {
  reject_experts = TRUE
  rejection_threshold = 0.05
  rejection_test = "distance_correlation"
  rej_expert_str <- ifelse(reject_experts, paste0("RejE(", rejection_threshold, ")"), "noR")
  output_file= decoupler_test_output_file_name(reject_experts, rejection_threshold, Sys.Date(), rel_dev = FALSE)

  min_nr_experts = 1
  studies <- load_data_49(relative_dev_folder=FALSE)
  run_decoupler_summarizer_dependence_analysis(
    studies,
    metric_list,
    output_file,
    reject_experts=reject_experts,
    rejection_threshold=rejection_threshold,
    rejection_test=rejection_test,
    min_nr_experts=min_nr_experts
  )
}
