% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pw_cum_mat.R
\name{pw_cum_mat}
\alias{pw_cum_mat}
\title{Generate Piecewise Constant Cumulative Hazard Matrix}
\usage{
pw_cum_mat(y, knots)
}
\arguments{
\item{y}{vector of event times.}

\item{knots}{increasing vector of cutpoints}
}
\value{
a numeric matrix, with rows corresponding to elements of y.
}
\description{
This helper function takes in a vector of event times, and
  generates a matrix dividing each time into a vector giving the time accrued within each
  interval between consecutive elements of a vector of knots.
}
