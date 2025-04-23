library(CopulaSEJ)
source("dev/dev_utils.R")

studies <- load_data_49(relative_dev_folder=FALSE) |> filter_studies_few_questions(min_questions = 11)
#studies <- studies |> filter_study_keep_ids(33)
studies_cross <- create_cross_validation_sets_list(studies)
summary_fun <- get_three_quantiles_summarizing_function()$f
cdf_decoupler <- get_CDF_decoupler()
hypothesis_tests <- list("classical_calibration", "distance_correlation")

# I would want to see the distribution of p values for the classical calibration score vs the distance correlation one.
study_ids <- c()
expert_ids <- c()
p_test <- c()
p_value <- c()
for (cross_study in studies_cross) {
  for (cross_i in nrow(cross_study)) {
    arr_format <- df_format_to_array_format(cross_study$training[[cross_i]],
                                            cross_study$test[[cross_i]],
                                            summary_fun)

    nr_experts <- arr_format$nr_experts
    study_ids <- append(study_ids, rep(arr_format$study_id, nr_experts * length(hypothesis_tests)))
    expert_ids <- append(expert_ids, rep(seq_len(nr_experts), length(hypothesis_tests)))
    for (test in hypothesis_tests) {
      p_test <- append(p_test, rep(test, nr_experts))
      test_res <- p_values_test(
        arr_format$training_summaries,
        arr_format$training_realizations,
        test,
        cdf_decoupler)
      p_value <- append(p_value, test_res)
    }
  }
}
df <- tibble::tibble(
  study_id = factor(study_ids, levels = unique(study_ids)),
  expert_id = expert_ids,
  test = p_test,
  p_value = p_value
)

df |> ggplot(aes(x=p_value, fill = test)) +
  geom_histogram(binwidth = 0.05, boundary=0) +
  theme_minimal() +
  labs(title = "Distribution of p-values for different tests",
       x = "p-value",
       y = "Count") +
  facet_wrap(~test, scales = "free_y") +
  scale_fill_manual(values=c("#FF9999", "#66B3FF"), guide="none")

