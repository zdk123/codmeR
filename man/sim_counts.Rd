\name{sim_counts}
\alias{sim_counts}
\title{Simulate count data from various theoretical distributions}
\usage{
sim_counts(n = 1, mu, Sigma = NULL, method = "ZINB", negscheme = exp,
  targetneff, ...)
}
\arguments{
  \item{n}{-> number of samples}

  \item{mu}{-> vector of 'OTU' means across samples}

  \item{Sigma}{-> positive definite matrix, one will be
  generated if not provided}

  \item{negscheme}{-> function to deal with negative values
  (for method MVLN data)}
}
\description{
Simulate count data from various theoretical distributions
}
\note{
required: length(mu) == nrow(Sigma) == ncol(Sigma)
}

