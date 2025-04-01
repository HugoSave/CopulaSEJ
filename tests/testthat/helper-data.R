get_scrambled_training_data <- function() {
  data.frame(
    `5th percentile` = c(70.84, 21.18, 15.00, 5.00, 65.00, 55.00),
    `50th percentile` = c(118.06, 35.30, 25.00, 10.00, 75.00, 65.00),
    `95th percentile` = c(177.10, 52.95, 27.00, 13.00, 77.00, 67.00),
    realization = c(213.78, -0.74, 12.00, 4.00, 49.00, 55.69),
    expert_id = c(2, 2, 1, 1, 1, 2),
    question_id = c(3, 1, 4, 3, 1, 4),
    study_id = c(1, 1, 1, 1, 1, 1),
    check.names=FALSE
  )
}

get_training_data <- function() {
  data.frame(
    `5th percentile` = c(70.84, 21.18, 15.00, 5.00, 65.00, 55.00),
    `50th percentile` = c(118.06, 35.30, 25.00, 10.00, 75.00, 65.00),
    `95th percentile` = c(177.10, 52.95, 27.00, 13.00, 77.00, 67.00),
    realization = c(213.78, -0.74, 12.00, 4.00, 49.00, 55.69),
    expert_id = c(1, 1, 1, 2, 2, 2),
    question_id = c(1, 3, 4, 1, 3, 4),
    study_id = c(1, 1, 1, 1, 1, 1),
    check.names=FALSE
  )
}

get_test_data <- function() {
  get_training_data() |> dplyr::filter(question_id == 3)
}

example_training_assessments_array <- function() {
  training_df <- get_training_data()
  data <- training_df |> dplyr::select(`5th percentile`, `50th percentile`, `95th percentile`)
  # split by quetsion
  data_matrix <- data |> split(training_df$question_id) |> abind::abind(along=3) |> aperm(c(3,1,2))
  dimnames(data_matrix) <- list(default_Q_names(training_df$question_id |> unique()), default_E_names(training_df$expert_id |> unique()), c("5th percentile", "50th percentile", "95th percentile"))
  data_matrix
}

example_test_matrix <- function() {
  get_test_data() |> dplyr::select(`5th percentile`, `50th percentile`, `95th percentile`) |> as.matrix()

}

example_error_distributions <- function(supports = NULL) {
  make_dist <- function(s1, s2, shape1=1.5, shape2=1.3) {
    force(s1)
    force(s2)
    force(shape1)
    force(shape2)
    list(
      pdf = function(x) extraDistr::dnsbeta(x, shape1=shape1, shape2=shape2, min=s1, max=s2),
      cdf = function(x) extraDistr::pnsbeta(x, shape1=shape1, shape2=shape2, min=s1, max=s2),
      support = c(s1,s2)
    )
  }
  if (is.null(supports)) {
    list(
      make_dist(0.5, 3),
      make_dist(0.8, 2),
      make_dist(0.2, 1.5),
      make_dist(0.1, 1),
      make_dist(0.1, 1),
      make_dist(0.1, 1)
    )
  } else {
    supports |> purrr::map(function(support) make_dist(support[1], support[2]))
  }
}

get_test_pseudo_obs <- function() {
  # these come from study 1 using ratio error metric and three quantiles
  pseudo_obs <- matrix(
    c(0.10, 0.1, 0.1, 0.9, 0.9, 0.1, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1,
      0.75, 0.8, 0.8, 0.5, 0.3, 0.4, 0.4, 0.7, 0.7, 0.3, 0.2, 0.4,
      0.75, 0.9, 0.9, 0.8, 0.8, 0.8, 0.6, 0.8, 0.9, 0.7, 0.8, 0.8,
      0.90, 0.7, 0.7, 0.7, 0.7, 0.5, 0.1, 0.4, 0.6, 0.6, 0.6, 0.7,
      0.50, 0.5, 0.5, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.2, 0.3, 0.2,
      0.60, 0.6, 0.6, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.5, 0.4, 0.3,
      0.20, 0.2, 0.2, 0.4, 0.5, 0.6, 0.5, 0.5, 0.5, 0.8, 0.7, 0.6,
      0.40, 0.4, 0.4, 0.3, 0.4, 0.9, 0.7, 0.9, 0.8, 0.9, 0.9, 0.9,
      0.30, 0.3, 0.3, 0.6, 0.6, 0.7, 0.8, 0.6, 0.4, 0.4, 0.5, 0.5),
    nrow = 9, byrow = TRUE
  )

  rownames(pseudo_obs ) <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
  pseudo_obs
}
