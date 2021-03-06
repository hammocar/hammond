% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{exactci}
\alias{exactci}
\title{Computes the Clopper/Pearon exact ci for a binomial success probability for x successes out of n trials with confidence coefficient conflev}
\usage{
exactci(x, n, conflev = 0.95)
}
\arguments{
\item{x}{Number of successes}

\item{n}{Number of trials}

\item{conflev}{Confidence level}
}
\value{
a list with lower bound, point estimate, and upper bound.
}
\description{
Computes the Clopper/Pearon exact ci for a binomial success probability for x successes out of n trials with confidence coefficient conflev
}
