% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_covariance.R
\name{make_ccov}
\alias{make_ccov}
\title{Compute the conditional covariances between effects in a nested model given observed data}
\usage{
make_ccov(data, group_effect = FALSE, ...)
}
\arguments{
\item{data}{An object of class \code{nesteddata}}

\item{group_effect}{A logical indicating wether the model includes an overall group random effect.}
}
\value{
A closure \code{ccov(effect1, effect2)} taking specifications of effects. The effects should
be specified as character vectors with the nesting structure of the effect in question. For example,

\code{ccov(c("sire1", "dam12"), "sire1")}

returns the conditional covariance matrix \code{Cov(beta_12, alpha_1 | y^o)}. A \code{NULL} value of either
effect corresponds to the overall group effect.
}
\description{
Compute the conditional covariances between effects in a nested model given observed data
}
