\name{sparseSigma}
\alias{sparseSigma}
\title{One way of generating a sparse, positive definite [square] correlation matrix}
\usage{
sparseSigma(dimension, sparsity = 0.6, center = 0)
}
\arguments{
  \item{dimension}{number of rows/cols}

  \item{sparsity}{controls probability of sparsity by
  setting the standard deviation of truncated normal
  distribution. Lower SD -> higher probability of a draw}

  \item{center}{correlations distributed around center
  [default = 0.0]}
}
\description{
One way of generating a sparse, positive definite [square]
correlation matrix
}

