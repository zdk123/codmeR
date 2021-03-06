\name{rmvpois}
\alias{rmvpois}
\title{Generate multivariate poisson data,
with counts approximately correlated according to Sigma}
\usage{
rmvpois(n, mu, Sigma, ...)
}
\arguments{
  \item{n}{number of samples to draw}

  \item{mu}{mean vector for variables (of length \eqn{D})}

  \item{Sigma}{\eqn{DxD} covariance or correlation matrix}

  \item{lambdas}{supply rate parameter (instead of mu)}

  \item{ps}{probability of zeros (instead of mu)}
}
\value{
\eqn{Dxn} matrix with zi-poisson data
}
\description{
Generate multivariate poisson data, with counts
approximately correlated according to Sigma
}

