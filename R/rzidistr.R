#' Draw samples from a zero-inflated poisson distribution
#'
#' @param n the number of samples to draw
#' @param lambda The poisson rate parameter
#' @param pstr0 probability of drawing a zero
#' @return Poisson counts of length \eqn{n}
#' @importFrom stats qpois dpois runif
#' @export
rzipois <- function(n, lambda, pstr0 = 0) {
    ans <- rpois(n, lambda)
    ans <- ifelse(runif(n) < pstr0, 0, ans)
    prob0 <- exp(-lambda)
    deflat.limit <- -1/expm1(lambda)
    ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
    if (any(ind0)) {
        pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * pstr0[ind0]
        ans[ind0] <- qpois(p = runif(sum(ind0), 
                    min = dpois(0, lambda[ind0])), lambda[ind0])
        ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
    }
    ans[pstr0 < deflat.limit] <- NaN
    ans[pstr0 > 1] <- NaN
    ans
}

#qzipois <- function(p, lambda, pstr0 = 0) {
#    LLL <- max(length(p), length(lambda), length(pstr0))
#    if (length(p) != LLL) 
#        p <- rep(p, len = LLL)
#    if (length(lambda) != LLL) 
#        lambda <- rep(lambda, len = LLL)
#    if (length(pstr0) != LLL) 
#        pstr0 <- rep(pstr0, len = LLL)
#    ans <- p
#    ans[p <= pstr0] <- 0
#    pindex <- (p > pstr0)
#    ans[pindex] <- qpois((p[pindex] - pstr0[pindex])/(1 - pstr0[pindex]), 
#        lambda = lambda[pindex])
#    deflat.limit <- -1/expm1(lambda)
#    ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
#    if (any(ind0)) {
#        pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * exp(-lambda[ind0])
#        ans[p[ind0] <= pobs0] <- 0
#        pindex <- (1:LLL)[ind0 & (p > pobs0)]
#        Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * exp(-lambda[pindex])
#        ans[pindex] <- qpospois((p[pindex] - Pobs0)/(1 - Pobs0), 
#            lambda = lambda[pindex])
#    }
#    ans[pstr0 < deflat.limit] <- NaN
#    ans[pstr0 > 1] <- NaN
#    ans[p < 0] <- NaN
#    ans[p > 1] <- NaN
#    ans
#}

#qpospois <- function (p, lambda) {
#    ans <- qpois(ppois(0, lambda, lower.tail = FALSE) * p + dpois(0, 
#        lambda), lambda = lambda)
#    ans[p > 1] <- NaN
#    ans[p < 0] <- NaN
#    ans[p == 1] <- Inf
#    ans
#}

#' @keywords internal
.zipois_getLam <- function(mu, S) {
    S <- max(sqrt(mu), S)
    (S^2/mu) + mu - 1
}

#' @keywords internal
.zipois_getP <- function(mu, S) {
    S <- max(sqrt(mu), S)
    (S^2 - mu) / (mu^2 - mu + S^2)
}

#' Generate multivariate, Zero-inflated poisson data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param lambdas supply rate parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzipois
#' @export
rmvzipois <- function(n, mu, Sigma, lambdas, ps, ...) {
    d   <- length(mu)
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))

    if (missing(lambdas) || missing(ps)) {
        lambdas <- unlist(lapply(1:length(SDs), function(i) .zipois_getLam(mu[i], SDs[i])))
        ps   <- unlist(lapply(1:length(SDs), function(i) .zipois_getP(mu[i], SDs[i])))
    }
    if (length(mu) == 1) stop("Need more than 1 variable")
    normd  <- rmvnorm(n, rep(0, d), Cor)
    unif   <- pnorm(normd)
    data <- t(VGAM::qzipois(t(unif), lambdas, pstr0=ps, ...))
    data <- .fixInf(data)
    return(data)
}

#rzinegbin <- function (n, size, prob = NULL, munb = NULL, pstr0 = 0)  {
#    if (length(munb)) {
#        if (length(prob)) 
#            stop("arguments 'prob' and 'munb' both specified")
#        prob <- size/(size + munb)
#    }
#    use.n <- if ((length.n <- length(n)) > 1) 
#        length.n
#    else if (!is.numeric(n, integer.valued = TRUE, allowable.length = 1, 
#        positive = TRUE)) 
#        stop("bad input for argument 'n'")
#    else n
#    pstr0 <- rep(pstr0, len = use.n)
#    size <- rep(size, len = use.n)
#    prob <- rep(prob, len = use.n)
#    ans <- rnbinom(n = use.n, size = size, prob = prob)
#    ans <- ifelse(runif(use.n) < pstr0, rep(0, use.n), ans)
#    prob0 <- rep(prob^size, len = use.n)
#    deflat.limit <- -prob0/(1 - prob0)
#    if (any(ind0, na.rm = TRUE)) {
#        pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
#        ans[ind0] <- rposnegbin(sum(ind0, na.rm = TRUE), size = size[ind0], 
#            prob = prob[ind0])
#        ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
#    }
#    ans[pstr0 < deflat.limit] <- NaN
#    ans[pstr0 > 1] <- NaN
#    ans
#}


#qzinegbin <- function (p, size, prob = NULL, munb = NULL, pstr0 = 0)  {
#    if (length(munb)) {
#        if (length(prob)) 
#            stop("arguments 'prob' and 'munb' both specified")
#        prob <- size/(size + munb)
#    }
#    LLL <- max(length(p), length(prob), length(pstr0), length(size))
#    if (length(p) != LLL) 
#        p <- rep(p, len = LLL)
#    if (length(pstr0) != LLL) 
#        pstr0 <- rep(pstr0, len = LLL)
#    if (length(prob) != LLL) 
#        prob <- rep(prob, len = LLL)
#    if (length(size) != LLL) 
#        size <- rep(size, len = LLL)
#    ans <- p
#    ind4 <- (p > pstr0)
#    ans[!ind4] <- 0
#    ans[ind4] <- qnbinom(p = (p[ind4] - pstr0[ind4])/(1 - pstr0[ind4]), 
#        size = size[ind4], prob = prob[ind4])
#    prob0 <- prob^size
#    deflat.limit <- -prob0/(1 - prob0)
#    ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
#    if (any(ind0)) {
#        pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
#        ans[p[ind0] <= pobs0] <- 0
#        pindex <- (1:LLL)[ind0 & (p > pobs0)]
#        Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
#        ans[pindex] <- qposnegbin((p[pindex] - Pobs0)/(1 - Pobs0), 
#            size = size[pindex], prob = prob[pindex])
#    }
#    ans[pstr0 < deflat.limit] <- NaN
#    ans[pstr0 > 1] <- NaN
#    ans
#}


#' @keywords internal
.zinegbin_getLam <- function(mu, S) {
    S   <- max(sqrt(mu)+1e-3, S)
    (mu + (mu^2 - mu + S^2) / mu) / 2
}

#' @keywords internal
.zinegbin_getP <- function(mu, lam) {
    1 - (mu / lam)
}

#' @keywords internal
.zinegbin_getK <- function(mu, S, lam) {
    S   <- max(sqrt(mu)+1e-3, S)
    (mu * lam) / (mu^2 - (mu * (lam + 1)) + S^2)
}


#' Generate multivariate, Zero-inflated negative binomial data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param munbs Rate/mean parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @param ks shape parameter
#' @param ... other arguments to the negative binomial distribution
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzinegbin
#' @export
rmvzinegbin <- function(n, mu, Sigma, munbs, ps, ks, ...) {
# Generate an NxD matrix of Zero-inflated poisson data,
# with counts approximately correlated according to Sigma
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    if (missing(munbs) || missing(ps) || missing(ks)) {
        if (missing(mu)) stop('mu is required')
        munbs <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getLam(mu[i], SDs[i])))
        ps   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getP(mu[i], munbs[i])))
        ks   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getK(mu[i], SDs[i], munbs[i])))
    }
    normd  <- rmvnorm(n, rep(0, length(munbs)), Sigma=Cor)
    unif   <- pnorm(normd)
    data <- t(VGAM::qzinegbin(t(unif), munb=munbs, size=ks, pstr0=ps, ...))
    data <- .fixInf(data)
    return(data)
}

#' @keywords internal
.fixInf <- function(data) {
    # hacky way of replacing infinite values with the col max + 1
    if (any(is.infinite(data))) {
       data <-  apply(data, 2, function(x) {
              if (any(is.infinite(x))) {
                   x[ind<-which(is.infinite(x))] <- NA
                   x[ind] <- max(x, na.rm=TRUE)+1
                 }
                x
                })
    }
    data
}


#' Draw samples from multivariate, correlated normal distribution
#' with counts correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param tol numerical tolerance for a zero eigenvalue (check for PD of Sigma)
#' @param empirical is Sigma the empirical correlation?
#' @return \eqn{Dxn} matrix with Gaussian data
#' @export
rmvnorm <- function(n=100, mu=rep(0,10), Sigma=diag(10), tol=1e-6, empirical=TRUE) {
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p))) 
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L]))) 
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0, nv = length(mu))$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    return(t(X))
}

#' Convert a symmetric correlation matrix to a covariance matrix
#' given the standard deviation
#'
#' @param cor a symmetric correlation matrix
#' @param sds standard deviations of the resulting covariance.
#' @return Covariance matrix of sample dimension as cor
cor2cov <- function(cor, sds) {
    if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
    cor * sds * rep(sds, each=nrow(cor))
}
