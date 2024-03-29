% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disclapmix_robust.R
\name{disclapmix_robust}
\alias{disclapmix_robust}
\title{Robust fitting}
\usage{
disclapmix_robust(
  x,
  clusters,
  rnd_obs = min(5L * clusters, round(nrow(x)/2)),
  rnd_tries = 10,
  ...
)
}
\arguments{
\item{x}{Dataset.}

\item{clusters}{The number of clusters/components to fit the model for.}

\item{rnd_obs}{Number of random observations in subset if initial fit fails}

\item{rnd_tries}{Number of tries with random subset if initial fit fails}

\item{\dots}{Passed on to `disclapmix()`}
}
\description{
A wrapper around `disclapmix()` that tries to avoid errors. 
Can sometimes avoid errors with SVD problems happening with 
`glm_method = 'internal_coef'` and `glm_method = 'internal_dev'`.
This is done by taking a random subset of observations.
}
\examples{
data(danes)
db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
fit <- disclapmix_robust(db, 3L)
fit

}
