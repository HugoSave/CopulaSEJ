% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copula_analysis.R
\name{calculate_marginal_posteriors}
\alias{calculate_marginal_posteriors}
\title{Calculates the unnormalized posteriors of $Q$ when conditioning on a single
$M_i^e$ random variable and assuming a flat prior.}
\usage{
calculate_marginal_posteriors(error_margins, m_matrix_observed, error_metric)
}
\arguments{
\item{error_margins}{d*E long list of error margins}

\item{m_matrix_observed}{Exd matrix of summary properties}

\item{error_metric}{Error metric of S3 class "error_metric"}
}
\description{
\deqn{
f_{Q\mid M_i^e = \dot{m_i^e}}(q) \propto f_{M_i^e \mid Q=q}(\dot{m_i^e} \mid q)
}
}
