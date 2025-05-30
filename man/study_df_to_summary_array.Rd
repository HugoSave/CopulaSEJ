% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_frame_format_utils.R
\name{study_df_to_summary_array}
\alias{study_df_to_summary_array}
\title{Converts a data frame to a nx3 matrix of assessments}
\usage{
study_df_to_summary_array(training_set, m_fun, k_percentiles = c(5, 50, 95))
}
\arguments{
\item{training_set}{data frame with columns expert_id, question_id, and the k_percentiles.}

\item{m_fun}{function}

\item{k_percentiles}{}
}
\value{
nxExd matrix
}
\description{
Converts a data frame to a nx3 matrix of assessments
}
