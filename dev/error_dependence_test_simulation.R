devtools::load_all("../")
library(dplyr)
library(energy)
source("dev_utils.R")


calculate_errors <- function(study_data, m_func, e_func, percentiles=c(5,50,95), flattened=TRUE) {
  if (nrow(study_data) == 0) {
    stop("Must have atleast one data row")
  }
  percentiles_col_names = k_percentiles_to_colname(percentiles)
  question_ids = sort(unique(study_data$question_id))
  expert_ids = sort(unique(study_data$expert_id))

  m_length = m_func(study_data[percentiles_col_names][1,]) |> length()
  nr_questions = length(question_ids)

  error_data = array(NA, dim = c(nr_questions, length(expert_ids), m_length))
  realizations = array(NA, dim = c(nr_questions, 1))

  for (q_id in question_ids){
    for (e_id in expert_ids){
      # get the data for the expert and question
      data = study_data |> filter(expert_id == e_id & question_id == q_id)
      if (nrow(data) != 1) {
        # print q_id and e_id
        stop(paste("Must have exactly one data row. Error for q:id", q_id, "e_id:", e_id))
      }

      assessments = data |> select(all_of(percentiles_col_names))
      m = m_func(assessments) # removes column names
      q = data$realization
      errors = as.vector(e_func(m, q), "double")
      error_data[q_id, e_id, ] = errors
      if (e_id == 1) {
        realizations[q_id, 1] = data$realization
      }
    }
  }
  if (flattened) {
    error_flattened <- matrix(aperm(error_data, c(1,3,2)), nrow=nr_questions)
    return(list(error_flattened=error_flattened, realizations=realizations))
  } else {
    return(list(error=error_data, realizations=realizations))
  }
}

error_difference <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  return(m - q)
}

error_difference_2 <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  return(q-m)
}

error_relative <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  return((q - m) / q)
}

error_relative_m_denom <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  return((q - m) / m)
}

error_symm_relative <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  ret_m = (q-m)/(m+q)
  # if some of the m elements are 0 while also q are 0, return 1 for those indices
  # if (q == 0) {
  #   m_zero_indices = which(m == 0)
  #   ret_m[m_zero_indices] = 1
  #   return(ret_m)
  # }
  return(ret_m)
}

# Log transform should not change anything
# error_log <- function(m, q){
#   # m is a D length vector. q is a scalar. Returns a D length vector
#   return(log(q/m))
# }
#
# error_log_2 <- function(m, q){
#   # m is a D length vector. q is a scalar. Returns a D length vector
#   return(log(m/q))
# }

error_ratio <- function(m, q){
  # m is a D length vector. q is a scalar. Returns a D length vector
  return(m/q)
}

# complete_data_df <- data_list_form |> do.call(what=rbind)

get_error_CDF_fun <- function(m_matrix, k_percentiles=c(5,50,95), overshoot=0.1){
  # fit distributions to m_matrix. m_matrix should be a matrix with E rows and D columns
  E = nrow(m_matrix)
  L = min(m_matrix)
  U = max(m_matrix)
  L_star = L - overshoot * (U - L)
  U_star = U + overshoot * (U - L)
  cdf_values <- c(0, k_percentiles/100, 1)
  # extended_m_matrix <- abind::abind(matrix(L_star, nrow=E, ncol=1), m_matrix, matrix(U_star, nrow=E, ncol=1), along=2)
  # distributions <- purrr::array_branch(extended_m_matrix, 1) |>
  #   purrr::map(\(fractiles) {
  #     linear_distribution_interpolation(fractiles, cdf_values)
  #   })

  CDF_error <- function(m, q){
    # m is a D length vector. q is a scalar. Returns a D length vector
    dist <- linear_distribution_interpolation(c(L_star, m, U_star), cdf_values)
    dist$cdf(q)
  }

  return(CDF_error)
}

run_error_dependence_simulation <- function(studies=NULL, metrics=NULL, output_file="output/dcors_df.rds") {
  if (is.null(studies)) {
    studies <- load_data_49()
  }
  if (is.null(metrics)) {
    metrics = list(
      list(m=m_three_quantiles, e=get_error_CDF_fun, name="3 quantiles and CDF error"),
      list(m=m_median, e=error_difference_2, name="median and q-m")
      #list(m=m_median, e=error_relative, name="median and (q-m)/q"),
      #list(m=m_median, e=error_relative_m_denom, name="median and (q-m)/m"),
      #list(m=m_median, e=error_symm_relative, name="median and (q-m)/(m+q)"),
      ##list(m=m_median, e=error_log, name="median and log(q/m)"),
      ##list(m=m_median, e=error_log_2, name="median and log(m/q)"),
      #list(m=m_median, e=error_ratio, name="median and m/q"),
      #list(m=m_three_quantiles, e=error_difference_2, name="3 quantiles and q-m"),
      ##list(m=m_three_quantiles, e=error_relative, name="3 quantiles and (q-m)/q"),
      ##list(m=m_three_quantiles, e=error_relative_m_denom, name="3 quantiles and (q-m)/m"),
      #list(m=m_three_quantiles, e=error_symm_relative, name="3 quantiles and (q-m)/(m+q)"),
      #list(m=m_three_quantiles, e=error_ratio, name="3 quantiles and m/q")
      ##list(m=m_3_quantiles, e=error_log, name="3 quantiles and log(q/m)")
    )
  }

  percentiles = c(5,50,95)
  dcors = list()
  for (study_id in seq_along(studies)) {
    print(paste("Running study:", study_id))
    combined_data <- studies[[study_id]] |> dplyr::bind_rows()
    # for each study_id and expert id, calculate the median. Output as matrix

    nr_questions = length(unique(combined_data$question_id))

    study_results = purrr::map(metrics, \(m_e) {
      if (m_e$name == "3 quantiles and CDF error") {
        e_fun <- m_e$e(combined_data)
        error_data = calculate_errors(combined_data, m_e$m, e_fun)
      } else {
        error_data = calculate_errors(combined_data, m_e$m, m_e$e)
      }
      # check if inf
      if (!all(is.finite(error_data$error_flattened))) {
        return(list(dcor=NA, name=m_e$name, study_id=study_id, dcor_test=NA, dcor_T_test=NA))
      }
      dcorT_result = energy::dcorT.test(error_data$error_flattened, error_data$realizations)
      dcor_result = energy::dcor.test(error_data$error_flattened, error_data$realizations, R=200)
      list(dcor=dcor_result$estimates["dCor"], name=m_e$name, study_id=study_id, dcor_test=dcor_result, dcor_T_test=dcorT_result)
    })
    # concat a and dcors
    dcors = c(dcors, study_results)
    # for question id
  }
  dcors_df = dcors |> purrr::map_dfr(~tibble(study_id=.x$study_id, dcor=.x$dcor, name=.x$name,
                                             dcor_test=list(.x$dcor_test),
                                             dcor_T_test=list(.x$dcor_T_test))) |> mutate(study_id=factor(study_id))
  saveRDS(dcors_df, output_file)
}




