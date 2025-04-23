# This files include utility functions that assume the data input is in a data frame format
# The required columns of the data frame is
#   study_id, question_id, expert_id
# together with percentile columns
# 'xth percentile' where x is an integer between 0 and 100


k_percentiles_to_colname <- function(k_percentiles) {
  # assert that we have a numeric vector with integers
  checkmate::assert_integerish(k_percentiles)
  paste0(k_percentiles, "th percentile")
}

#' Converts a data frame to a list of arrays suitable for the CopulaSEJ package
#'
#'
#' @param training_set,test_set data frames with columns expert_id, question_id, and the k_percentiles.
#' @param m_fun function to apply to the k_percentiles values to get summary values
#' @param percentiles a vector of percentiles 'k' to be used in the data frame.
#'
#' @returns a list with three elements:
#' - training_summaries: a 3D array of shape (n, E, d) with the training set summaries
#' - training_realizations: a vector of length n with the training set realizations
#' - test_summaries: a 2D array of shape (E, d) with the test set summaries
#' @export
#'
#' @examples
df_format_to_array_format <- function(training_set, test_set, m_fun, percentiles=c(5,50,95)) {
  checkmate::assert_data_frame(training_set, any.missing = FALSE)
  checkmate::assert_data_frame(test_set, any.missing = FALSE)
  checkmate::assert_subset(c("question_id", "expert_id"), colnames(training_set))
  checkmate::assert_subset(c("question_id", "expert_id"), colnames(test_set))

  training_summaries <- study_df_to_summary_array(training_set, m_fun, percentiles)
  training_realizations <- study_df_to_realizations(training_set)
  test_summaries <- study_df_single_question_to_summary_matrix(test_set, m_fun, percentiles)
  checkmate::assert_set_equal(dim(training_summaries)[2:3], dim(test_summaries))

  ret_list <- list(
    training_summaries = training_summaries,
    training_realizations = training_realizations,
    test_summaries = test_summaries,
    nr_experts = dim(training_summaries)[2],
    nr_questions = dim(training_summaries)[1],
    d = dim(training_summaries)[3]
  )

  # check if study_id is in the data frame
  study_id <- get_study_id_if_exists(training_set)
  if (!is.null(study_id)) {
    # add study_id to the list
    ret_list$study_id <- study_id
  }

  ret_list
}

get_study_id_if_exists <- function(df) {
  # check if study_id is in the data frame
  selection <- df |> dplyr::select(dplyr::any_of("study_id")) |> dplyr::distinct()
  stopifnot(nrow(selection) <= 1)
  if (nrow(selection) == 1 ) {
    selection$study_id[[1]]
  } else {
    NULL
  }
}

#' Extracts the percentile columns from the dataframe and returns a matrix
#'
#' @param study_data - a data frame with percentile columns. n rows
#' @param percentiles - the percentiles present in the data frame. Length d. Integer
#'
#' @returns - a matrix (nxd) with the assessment data
#' @export
#'
percentiles_from_dataframe <- function(df, percentiles) {
  col_names <- k_percentiles_to_colname(percentiles)
  assertthat::assert_that(all(col_names %in% colnames(df)))
  row_names<-rownames(df)

  as.matrix(df[col_names], dimnames=list(row_names, col_names))
}

create_cross_validation_sets_list <- function(study_list) {
  study_list |> purrr::map(create_cross_validation_sets)
}

create_cross_validation_sets <- function(study_data) {
  # Nests by question id and then creates all possible combinations of training and test sets
  fold_combinations <- study_data |> tidyr::nest(.by = "question_id",
                                                 .key = "question_data") |>
    ExhaustiveFolds(1) |> dplyr::mutate(
      training = purrr::map(training, unnest, "question_data"),
      test = purrr::map(test, unnest, "question_data")
    )
  # returns a tibble with columns question_id, training, test. training and test are themselves also tibbles.
  fold_combinations
}

#' Converts a data frame to a nx3 matrix of assessments
#'
#' @param training_set data frame with columns expert_id, question_id, and the k_percentiles.
#' @param m_fun function
#' @param k_percentiles
#'
#' @returns nxExd matrix
#' @export
#'
study_df_to_summary_array <- function(training_set, m_fun, k_percentiles=c(5,50,95)) {
  # assert only one study
  col_names <- colnames(training_set)
  if ("study_id" %in% col_names) {
    assertthat::are_equal(length(unique(training_set$study_id)), 1)
  }
  assertthat::assert_that("question_id" %in% colnames(training_set), "expert_id" %in% colnames(training_set))
  # Convert the training set to a matrix
  ordered_training <- training_set |> dplyr::arrange(question_id, expert_id)
  assessments <- percentiles_from_dataframe(ordered_training, k_percentiles)
  if (is.null(colnames(assessments))) {
    colnames(assessments) <-  k_percentiles %>% paste0("M", .)
  }
  m <- m_fun(assessments) # dim = (n*E, d)

  m_3d <- m |> split.data.frame(ordered_training$question_id) |>
    abind::abind(along = 3) |> # dim = (E, d, n)
    aperm(c(3, 1, 2)) # dim = (n, E, d)

  dimnames_1 <- ordered_training$question_id |> unique()  %>% paste0("Q", .)
  dimnames_2 <- ordered_training$expert_id |> unique() %>% paste0("E", .)
  dimnames(m_3d)[1:2] <- list(dimnames_1, dimnames_2)

  m_3d
}


study_df_to_realizations <- function(training_set) {
  # assert only one study
  col_names <- colnames(training_set)
  if ("study_id" %in% col_names) {
    assertthat::are_equal(length(unique(training_set$study_id)), 1)
  }
  checkmate::expect_subset(c("question_id", "realization"), colnames(training_set))
  # Convert the training set to a matrix
  summarised <- training_set |> dplyr::group_by(question_id) |>
    dplyr::summarise(realization = realization |> unique(),
                     question_id = unique(question_id))
  summarised |> dplyr::arrange(question_id) |> dplyr::pull(realization)
}

study_df_single_question_to_summary_matrix <- function(test_set, m_fun, k_percentiles) {
  # Convert the test set to a matrix
  arr <- study_df_to_summary_array(test_set, m_fun, k_percentiles)
  checkmate::assert_true(dim(arr)[1] == 1)
  abind::adrop(arr, 1) # if m_fun is median then the third dimension is also 1 so in that case we need to make sure not to drop that
}
