#' compute p-vals of proportionality (two ways) via aitchison variance
#'
#' @param data OTUs/features in columns, samples are rows
#' @export
proportionality_filter <- function(data, normfun=function(x) x, nboots=1333,
                                    mar=2, n.core=1, ...) {
    
    .astat <- function(data, indices)   triu(aitvar(data[indices,,drop=FALSE]))
    .pstat <- function(data, indices)   triu(aitvar(apply(data[indices,], 2, 
                                                function(x) sample(x))))

    if (mar == 1) data <- t(data)
    d <- ncol(data)
    n <- nrow(data)
    pobjs <- tavboot(data, .astat, .pstat, R=nboots, ...)
    pvals <- pval(pobjs)
    pmat    <- triu2diag(pvals)
    rownames(pmat) <- colnames(pmat) <- colnames(data)
    pmat
}


triu <- function(x) x[upper.tri(x)]
tril <- function(x) x[lower.tri(x)]

triu2diag <- function(x, diagval=0) {
    e <- length(x)
    n <- .5 * (sqrt(8*e + 1)+1)
    mat <- matrix(0, n, n)
    mat[lower.tri(mat)] <- x
    mat <- mat + t(mat)
    diag(mat) <- diagval
    mat
}

pval <- function(x) {
    UseMethod('pval')
}

#' @keywords internal
#' @importFrom boot boot
tavboot <- function(data, statisticboot=astat, statisticperm, R, ...) {
    res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=3, ...)
    null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=3)
    class(res) <- 'list'
    structure(c(res, list(null_av=null_av)), class='tavboot')
}

#' @keywords internal
pval.tavboot <- function(x, sided='both', mar=2) {
# calculate 1 or 2 way pseudo p-val from boot object
# Args: a boot object
    if (sided != "both") stop("only two-sided currently supported")
    nparams  <- ncol(x$t)
    tmeans   <- colMeans(x$null_av$t);
    # calc whether center of mass is above or below the mean
    bs_above <- unlist(lapply(1:nparams, function(i) 
                    length(which(x$t[, i] > tmeans[i]))))
    pval <- ifelse(bs_above > x$R/2, 2*(1-bs_above/x$R), 2*bs_above/x$R)
    gt1  <- which(pval > 1)
    if (length(gt1) > 0) pval[gt1] <- 1
   # triu2diag(pval)
   pval
}

#' @export
aitvar <- function(x, ...) {
    UseMethod('aitvar')
}

#' @export
aitvar.default <- function(x, y) {
    if (length(x) != length(y)) stop('Error: data must be of same dimension')
    var(log(x/y))
}

#' compute aitchison variation of a matrix
#'
#' @import Rcpp
#' @export
aitvar.matrix <- function(x, mar=2) {
    if (!("fastaitvar" %in% ls())) sourceCpp("fastaitvar.cpp")
    avmat <- fastaitvar(x)
    colnames(avmat) <- rownames(avmat) <- colnames(x)
    avmat
}

#' @export
aitvar.data.frame <- function(x, mar=2, direction='pos') {
    aitvar(as.matrix(x))
}



