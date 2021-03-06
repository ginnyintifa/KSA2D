% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EB_diff_1Ddiff_wrap_functions.R
\name{comparison_time_points_1d}
\alias{comparison_time_points_1d}
\title{Compare between 2 conditions (paired) in 1D}
\usage{
comparison_time_points_1d(
  s1_col_name,
  s2_col_name,
  d1_data,
  nna_cutoff,
  all_seeds,
  permute_times,
  working_dir,
  compare_name
)
}
\arguments{
\item{s1_col_name}{names of the columns in the first condition}

\item{s2_col_name}{names of the columns in the second condition}

\item{d1_data}{dataframe for the data, rows are genes, columns are patients in 2 conditions}

\item{nna_cutoff}{genes with values in at least the number of samples are included in comparison}

\item{all_seeds}{seeds for the permutation}

\item{permute_times}{number of permutations}

\item{working_dir}{the directory output files will be deposited in}

\item{compare_name}{a label for the anlaysis}
}
\description{
Compare between 2 conditions (paired) in 1D
}
\examples{
camparison_time_points_1d
}
\keyword{1D}
