#' Draw mean vector from a log-normal vector & find the logsd parameter that
#' matches the desired n_eff
#'
#' @export
log_norm_mean_vector <- function(D, logmean=1, neff) {
# Generate an micrbial community-like vector of D means
    if (neff > D) stop("Effective species number is higher than number of species")
    logsd <- optim(0, fitfn_lnorm, neff=neff, D=D, logmean=logmean, method="Brent", lower=0, upper=1e1, control = list(abstol=1e-8))
    means <- replicate(100, rlnorm(D, logmean, logsd$par))
    neffs <- apply(means, 2, neff)
    ind   <- which.min(abs(neff - neffs))
    means[,ind]
}

#' @keywords internal
fitfn_lnorm <- function(logsd, neff, logmean, D) {
    abs(neff - mean(replicate(10, neff(rlnorm(D, logmean, logsd)))))
}

#' Draw mean vector from a vector of fixed means & pick the abundance values for a
#' fixed number of bugs with higher means that also matches the desired n_eff
#'
#' @export
fixed_mean_vector <- function(D, neff, val=2, numhigher=ceiling(D/50)) {
    if (neff > D) stop("Effective species number is higher than number of species")
    y <- optim(rep(val, numhigher), fitfn_fixed, x=(x<-rep(val, D-numhigher)), neff=neff, method="L-BFGS-B", lower=val, upper=Inf, control = list(factr = .5))
    sample(c(x,y$par))
}

#' @keywords internal
fitfn_fixed <- function(y, x, neff) {
    abs(neff - neff(c(y,x)))
}


#' Simulate count data 
#'
#'   @param n -> number of samples
#'   @param mu -> vector of 'OTU' means across samples
#'   @param Sigma -> positive definite matrix, one will be generated if not provided
#'   @param negscheme -> function to deal with negative values
#'   @note required: length(mu) == nrow(Sigma) == ncol(Sigma)
#'   @export
sim_counts <- function(n = 1, mu, Sigma=NULL, method='ZIP', negscheme=exp, targetneff, ...) {

if (method == 'MVLN') { 
    require(MASS)
    corpar=as.matrix(Sigma)
    Y <- negscheme(mvrnorm(n, mu, Sigma))

} else if (method == "Poi") {

    require(corcounts)
    margins <- rep("Poi", length(mu))
    corstr <- 'unstr'
    
    if (length(Sigma) == 0) {
        corpar <- unstructured(length(mu))
    } else {
        corpar <- as.matrix(Sigma)
    }

    Y <- rcounts(N=n, margins=margins, mu=mu, corstr=corstr, corpar=corpar)
    
} else if (method == "MN") {
    Y <- negscheme(multivnomial(n, mu, as.matrix(Sigma)))
    
} else if (method == "ZIP") {
    Y <- zi_distr("zipois", n, mu, Sigma, ...)
    
} else if (method == "ZINB") {
    Y <- zi_distr("zinegbin", n, mu, Sigma, ...)
} else {
    stop("Error: method not supported")
}

    list(data=as.matrix(Y), empcor=cor(Y), incor=as.matrix(Sigma))
}


multivnomial <- function(n, mu, Sigma) {
    N <- sum(mu)
    data <- rmultinom(n, N, mu/N)
    c <- as.matrix(chol(Sigma))
    t(data) %*% c
}

#' @export
unsparseSigma <- function(dimension, avecor=0.5, sparsity=0.5) {

    require(msm) ; require(Matrix)
    cormat <- matrix(0, dimension, dimension)
    d <- dimension

    cors <- rtnorm((d*(d-1))/2, mean=avecor, lower=-.999, upper=.999, sd=sparsity)
    cormat[lower.tri(cormat, diag=F)] <- cormat[upper.tri(cormat, diag=F)] <- cors
    diag(cormat) <- 1
    return(nearPD(cormat, corr=TRUE)$mat)
}

#' @keywords internal
addmin <- function(x, jitter=1e-3) {

    x + abs(min(x)) + jitter

}

#' compute the shannon entropy from a vector (normalized internally)
#'
#' Shannon entropy is:
#'     sum [ x_i log(x_i) ]
#'      
#' @param x data vector
#' @return shannon entropy in base e
#' @export
shannon <- function(x) {
    x.f <- (x+1)/sum(x+1)
    -sum(x.f*log(x.f))
}

#' N_effective: Compute the exponential of the shannon entropy. linearizes shannon entropy, for a better diveristy metric (effective number of species)
#'
#' @param x data vector
#' @return N_eff in base e
#' @export
neff <- function(x) exp(shannon(x))

#' One way of generating a sparse, positive definite [square] correlation matrix
#'
#' @param dimension number of rows/cols
#' @param  sparsity controls probability of sparsity by setting the standard deviation of truncated normal distribution. Lower SD -> higher probability of a draw
#' @param center correlations distributed around center [default = 0.0]
#' @export
sparseSigma <- function(dimension, sparsity=0.6, center=0.0) {

    require(msm) ; require(Matrix)
    temp <- matrix(0, dimension, dimension)
    diag(temp) <- 1

    prob <- dimension / dimension^2

    for (i in 2:dimension) {
        for (j in 1:(i - 1)) {
            if ( sample(c(0,1), 1, prob=c(prob, 1-prob)) == 0 ) {
                temp[i, j] <- rtnorm(1, mean=center, lower=-.99, upper=.99, sd=sparsity)
                temp[j, i] <- temp[i, j]
            }
        }
    }
    return(as.matrix(nearPD(temp, corr=TRUE)$mat))
}


#' Jonathan Friedman's described method of computing a correlation matrx
#'
#' @param D Dimension
#' @param p probability that any otu pairs are perfectly correlated
#' @param var individual variance (same for all OTUs)
#' @param center average correlation (default 0.0)
sparseSigma_JF <- function(D, p, var=1, center=0.0) {
    require(Matrix)
    Sigma <- diag(D) + center
    diag(Sigma) <- diag(Sigma) - center
    n <- (D*(D-1))/2
    Sigma[lower.tri(Sigma)] <- sample(c(0,1), n, prob=c(1-p,p), replace=TRUE)
    signs <- sample(c(-1,1), n, prob=c(.5,.5), replace=TRUE)
    Sigma[lower.tri(Sigma)] <- Sigma[lower.tri(Sigma)] * signs
    Sigma[upper.tri(Sigma)] <- Sigma[lower.tri(Sigma)]
    return(nearPD(Sigma, corr=TRUE)$mat)
}


#' keywords internal
convert_MD_format <- function(matrix, textname="otus") {
    if (is.null(colnames(matrix))) colnames(matrix) <- 1:ncol(matrix)
    if (is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
    header <- c(textname, colnames(matrix))
    newmat <- cbind(rownames(matrix), matrix)
    newmat <- rbind(header, newmat)
    return(newmat)
}

