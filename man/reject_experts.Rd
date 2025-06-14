% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expert_rejection.R
\name{reject_experts}
\alias{reject_experts}
\title{Title}
\usage{
reject_experts(
  training_assessments,
  training_realizations,
  rejection_level = 0.1,
  test = "kruskal",
  decoupler = NULL,
  min_nr_experts = NULL
)
}
\arguments{
\item{training_realizations}{Q long vector with the realizations of the training data.}

\item{rejection_level}{numeric value between 0 and 1. The level of rejection.}

\item{test}{string indicating the test to use.}

\item{decoupler}{decoupler object. Only used if test is "kruskal".}

\item{min_nr_experts}{minimum number of experts to accept. If the number of
accepted experts is less than this value, the best experts are selected based on p-values.}

\item{training_estimates}{QxExD array with the training data.}
}
\value{
A list with the following elements:
\itemize{
\item accepted_estimates: a Qx\tilde{E}xD array with the accepted estimates
\item p_values: a vector with the p-values for each expert
}
}
\description{
Title
}
