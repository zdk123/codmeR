#' Sparse/penalized estimators of covariance matrices
#'
#' @importFrom huge huge huge.npn
#' @importFrom flare flare.tiger
#' @export
sparseiCov <- function(data, method, npn=FALSE, verbose=FALSE, ...) {

    if (npn) data <- huge.npn(data, verbose=verbose)

    args <- list(...)
    if (is.null(args$lambda)) {
        # TODO: check if these lambdas are correct for flare methods
        S <- cor(data)
        lamMax <- max(max(abs(triu(S)))+1e-5, 1e-5)
        lamMin <- max(min(abs(triu(S)))/4, 1e-5)

        lambda <- exp(seq(log(lamMax), log(lamMin), 
                    length.out = ifelse(is.null(args$nlambda), 10, args$nlambda)
                      ))
        args$lambda <- lambda
    }
    method <- switch(method, glasso = "glasso", mb = "mb",
                      tiger = "slasso", clime = "clime", ADM = "ADM",
                      stop("Method not supported"))

    if (method %in% c("glasso", "mb")) {
         do.call(huge::huge, c(args, list(x=data, method=method, verbose=verbose, 
                     cov.output = TRUE)))

    } else if (method %in% c("slasso", "clime")) {
         do.call(flare::flare.tiger, c(args, list(data=data, method=method,  
                     verbose=verbose)))
    
    } else if (method == "ADM") stop("support coming soon")
}
