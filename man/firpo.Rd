% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.r
\name{firpo}
\alias{firpo}
\title{Firpo estimator (Firpo)}
\usage{
firpo(y, t, Q, g, q)
}
\arguments{
\item{y}{Vector with outcome values.}

\item{t}{Vector with binary treatment indicator.}

\item{Q}{Conditional outcome distribution estimate. This should come in the for of an n x p matrix, where each column represents a conditional quantile. (See Kang & Schafer example)}

\item{g}{Propensity score.}

\item{q}{Quantile to be computed (e.g., q = 0.5 for the median.)}
}
\value{
A point estimate
}
\description{
Firpo estimator (Firpo)
}
