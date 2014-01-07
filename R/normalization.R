#####################################################################
# Different normalization schemes for microbiome counts (real or fake)
#
# @author Zachary Kurtz
# @date 10/10/2013
#####################################################################


#' @export
norm_pseudo  <- function(x) norm_to_total(x+1)

#' @export
norm_diric   <- function(x, rep=1) {
    require(VGAM)
    dmat <- rdiric(rep, x+1)
    norm_to_total(colMeans(dmat))
}

#' @export
norm_to_total <- function(x) x/sum(x)

#' @export
clr <- function(x, ...) {
    UseMethod('clr')
}

#' @keywords internal
clr.default <- function(x.f, tol=.Machine$double.eps) {
# Center-log-ratio transform
# Args:
#   x.f -> component vector
#   tol -> definition of machine zero
    zero <- (x.f >= tol)
    LOG <- log(ifelse(zero, x.f, 1))
    ifelse(zero, LOG - mean(LOG)/mean(zero), 0.0)
}

#' @keywords internal
clr.matrix <- function(x.f, mar=2, ...) {
    apply(x.f, mar, clr, ...)
}


#x_list <- vector('list', 1000)
#for (l in 1:1000) {
#    x_list[[l]] <- norm_diric(x, rep=l)
#}

#gmean <- do.call(function(x) exp(mean(log(x))), x_list)

#' @keywords internal
CSS <- function(x, ...) {
    UseMethod("CSS")
}
#' @keywords internal
CSS.default <- function(x, p=0.05, sl=1000) {
# Cumulative sum scaling Normalization Paulson et al 2013 (Nature Methods)
    xx <- x
    xx[x==0] <- NA
    qs <- quantile(xx, p=p, na.rm=TRUE)
    xx <- x - .Machine$double.eps
    normFactor <- sum(xx[xx <= qs])
    (x/normFactor)*sl
}
#' @keywords internal
CSS.matrix <- function(x, p=CSSstat(x), sl=1000, mar=2) {
    apply(x, mar, CSS, p=p, sl=sl)
}
#' @keywords internal
CSSstat <- function(mat, rel=0.1) {
    smat <- sapply(1:ncol(mat), function(i) {
        sort(mat[, i], decreasing = FALSE)
    })
    ref <- rowMeans(smat)
    yy  <- mat
    yy[yy == 0] <- NA
    ncols <- ncol(mat)
    refS  <- sort(ref)
    k     <- which(refS > 0)[1]
    lo    <- (length(refS) - k + 1)
    diffr <- sapply(1:ncols, function(i) {
            refS[k:length(refS)] - quantile(yy[, i], p = seq(0, 
                1, length.out = lo), na.rm = TRUE)})
    diffr2 <- apply(abs(diffr), 1, median, na.rm = TRUE)
    x <- which(abs(diff(diffr2))/diffr2[-1] > rel)[1]/length(diffr2)
    names(x) <- NULL
    x
}

#' @keywords internal
DESeq <- function(x, ...) {
    UseMethod("DESeq")
}
#' @keywords internal
DESeq.matrix <- function(mat, c) {
    # compute geometric mean along columns
    matt  <- mat
    matt[matt == 0] <- NA
    k_ref <- apply(matt, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
    krefmat <- matrix(rep(k_ref,ncol(mat)), nrow=nrow(mat))
    s_hat   <- apply(matt/krefmat, 2, median, na.rm=TRUE)
    if (missing(c)) {
        fn <- function(c, s_hat) abs(sum(log(c*s_hat)))
        c  <- optimize(fn, interval=0:10, s_hat=s_hat)$minimum
    }
    s <- c * s_hat
    smat <- matrix(rep(s,nrow(mat)), ncol=ncol(mat), byrow=TRUE)
    mat/smat
}



