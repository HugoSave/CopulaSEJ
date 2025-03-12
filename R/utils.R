

sample_log_unnormalized_density <- function(log_density, support, num_samples) {
  #MCMCpack::MCMCmetrop1R(log_density, (support[1] + support[2])/2, mcmc=num_samples)$batch
  #as.vector(mcmc::metrop(log_density, (support[1] + support[2])/2, num_samples)$batch)
  armspp::arms(num_samples, log_density, support[1] , support[2])
}

flatten_matrix_row_by_row <- function(matrix) {
  as.vector(t(matrix))
}


#' Splits data frame by questions and structures error order
#'
#' Takes a data frame with the columns expert_id, question_id, and the k_percentiles
#' and splits the data frame by question_id (in ascending order). For each sub
#' data frame the errors are calculated and put into a d*E long vector ordered
#' such that the errors from a single expert are kept adjacent. These vectors are then
#' put row by row in a matrix resulting in a nx(d*E) matrix where n is the
#' number of questions.
#'
#' @param dataframe frame data frame with columns expert_id, question_id, and the k_percentiles
#' @param error_metric error metric object
#' @param summarizing_function summarizing function
#' @param k_percentiles vector of percentiles
#'
#' @returns List of vectors
#'
split_dataframe_to_error_observations <- function(dataframe,
                                                     error_metric,
                                                     summarizing_function,
                                                     k_percentiles = c(5, 50, 95)) {
  # Split the training set into error observations
  # training set contains n*E rows where n is the number of questions and E is the number of experts
  ordered_training <- dataframe |> dplyr::arrange(question_id, expert_id)
  assessments <- data_frame_to_assessment_matrix(ordered_training, k_percentiles)
  m <- summarizing_function(assessments) # dim = (n*E, d)
  errors <- error_metric$f(m, ordered_training$realization) # dim = (n*E, d)
  # orders the second dimension such that a the errors for a single error type are kept adjacent
  # that is, after splitting by question id, the smaller error matrix is flattened in row order
  flattened_errors <- errors |> split.data.frame(ordered_training$question_id) |>
    purrr::map(flatten_matrix_row_by_row) |> do.call(what = rbind) # dim = (n, d*E)
  flattened_errors
}
