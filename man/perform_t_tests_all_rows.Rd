% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edec_stage_0.R
\name{perform_t_tests_all_rows}
\alias{perform_t_tests_all_rows}
\title{T-tests for comparing all rows of two matrices}
\usage{
perform_t_tests_all_rows(data_group_1, data_group_2)
}
\arguments{
\item{data_group_1}{A matrix of numbers.}

\item{data_group_2}{A matrix of numbers.}
}
\value{
A matrix where each row corresponds to comparisons between the
  corresponding rows in \code{data_group_1} and \code{data_group_2}. The
  first column contains the p-values, the second column contains the means of
  data_group_1 for each row, and the third column contains the means of
  data_group_2 for each row.
}
\description{
\code{perform_t_tests_all_rows} performs two sample T-tests comparing the
values in each row of one matrix against the values of the matching row in a
second matrix.
}

