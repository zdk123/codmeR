#' Sparse/penalized estimators of covariance matrices
#'
#' @importFrom huge huge huge.npn
#' @importFrom flare flare.tiger
#' @export
sparseiCov <- function(data, method, npn=FALSE, verbose=FALSE, ...) {

    if (npn) data <- huge.npn(data, verbose=verbose)
#    Cov    <- cov(data)
    method <- switch(method, glasso = "glasso", mb = "mb", ct = "ct",
                      tiger = "slasso", clime = "clime", ADM = "ADM",
                      stop("Method not supported"))

    if (method %in% c("glasso", "mb", "ct")) {
        huge::huge(data, method=method, verbose=verbose, 
                    cov.output = TRUE, ...)
    } else if (method %in% c("slasso", "clime")) {
        flare::flare.tiger(data, method=method, verbose=verbose, ...)
    } else if (method == "ADM") stop("support coming soon")
}
