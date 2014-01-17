#' compute p-vals of proportionality (two ways) via aitchison variance
#'
#' @param data OTUs/features in columns, samples are rows
#' @export
proportionality_filter <- function(data, normfun=function(x) x, nboots=1333,
                                    mar=2, n.core=1, ...) {
    
    .astat <- function(data, indices)   triu(aitvar(data[indices,,drop=FALSE]))
    .pstat <- function(data, indices)   triu(aitvar(apply(data[indices,], 2, function(x) sample(x))))

    if (mar == 1) data <- t(data)
    d <- ncol(data)
    n <- nrow(data)
    pobjs <- tavboot(data, .astat, .pstat, R=nboots, ncpus=n.core, ...)
    pvals <- pval(pobjs)
    pmat    <- triu2diag(pvals$pvals)
    rownames(pmat) <- colnames(pmat) <- colnames(data)
    structure(list(T.signed=triu2diag(pvals$sT), pvals=pmat), class='TAV')
}


triu <- function(x) x[upper.tri(x)]
tril <- function(x) x[lower.tri(x)]

triu2diag <- function(x, diagval=0) {
    e <- length(x)
    n <- .5 * (sqrt(8*e + 1)+1)
    mat <- matrix(0, n, n)
    mat[upper.tri(mat)] <- x
    mat <- mat + t(mat)
    diag(mat) <- diagval
    mat
}

#' @export
pval <- function(x) {
    UseMethod('pval')
}

#' @keywords internal
#' @importFrom boot boot
tavboot <- function(data, statisticboot, statisticperm, R, ncpus=1, ...) {
    res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=ncpus, ...)
    null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus)
    class(res) <- 'list'
    structure(c(res, list(null_av=null_av)), class='tavboot')
}

#' @keywords internal
pval.tavboot <- function(x, sided='both', mar=2) {
# calculate 1 or 2 way pseudo p-val from boot object
# Args: a boot object
    if (sided != "both") stop("only two-sided currently supported")
    nparams  <- ncol(x$t)
    tmeans   <- colMeans(x$null_av$t)
#    check to see whether Aitchison variance is unstable -- confirm 
#    that sample Aitchison variance is in 95% confidence interval of 
#    bootstrapped samples
    niters   <- nrow(x$t)
    ind95    <- max(1,round(.025*niters)):round(.975*niters)
    boot_ord <- apply(x$t, 2, sort)
    boot_ord95 <- boot_ord[ind95,]
    outofrange <- unlist(lapply(1:length(x$t0), function(i) {
            aitvar <- x$t0[i]
            range  <- range(boot_ord95[,i])
            range[1] > aitvar || range[2] < aitvar
        }))
    # calc whether center of mass is above or below the mean
    bs_above <- unlist(lapply(1:nparams, function(i) 
                    length(which(x$t[, i] > tmeans[i]))))
    is_above <- bs_above > x$R/2
    signedAV <- x$t0
    signedAV[is_above] <- -signedAV[is_above]
    pvals    <- ifelse(is_above, 2*(1-bs_above/x$R), 2*bs_above/x$R)
    pvals[pvals > 1]  <- 1
    pvals[outofrange] <- NaN
    list(sT=signedAV, pvals=pvals)
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
    if (mar == 2) x <- t(x)
    avmat <- fastaitvar(x)
    colnames(avmat) <- rownames(avmat) <- rownames(x)
    avmat
}

#' @export
aitvar.data.frame <- function(x, mar=2, direction='pos') {
    aitvar(as.matrix(x), mar)
}

#' @export
hamming.dist <- function(x, ...) {
    UseMethod('hamming.dist')
}

#' @export
hamming.dist.default <- function(x, y) {
    sum(x != y)
}

#' @export
hamming.dist.matrix <- function(x, mar=2) {
    if (mar == 1) x <- t(x)
    HamdistMat(x)
}

#' @export
hamming.dist.data.frame <- function(x, mar=2) {
    hamming.dist(as.matrix(x), mar)
}

#' @export
cooccurance <- function(x, ...) {
    UseMethod('cooccurance')
}

#' @export
cooccurance.default <- function(x, y) {
    sum((sign(x) + sign(y)) == 2)
}

#' @export
cooccurance.matrix <- function(x, mar=2) {
    if (mar == 1) x <- t(x)
    cooccurMat(x)
}


#' @export
coabsence <- function(x, ...) {
    UseMethod('coabsence')
}

#' @export
coabsence.default <- function(x, y) {
    sum((sign(x) + sign(y)) == 0)
}

#' @export
coabsence.matrix <- function(x, mar=2) {
    if (mar == 1) x <- t(x)
    coabsenceMat(x)
}

