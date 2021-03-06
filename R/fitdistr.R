#' from count data (ex HMP) fit parameters to OTU margins
#' and simulate a new community with those properties
#'
#' @param comm community: matrix of counts
#' @param distr distribution to fit (see fitdistr)
#' @param Sigma covariance structure (defaults to empirical cov of comm)
#' @param params optionally supply already fitted parameters
#' @param n number of samples (defaults to comm samples)
#' @param retParams if TRUE, return the fitted parameters
#' @param ... additional parameters to parameter fitting
#' @return community 
#' @export
synth_comm_from_counts <- function(comm, distr, Sigma=cov(t(comm)), params, n=ncol(comm), retParams=FALSE, ...) {
    D <- nrow(comm)
    if (missing(params))  params <- get_comm_params(comm, distr, ...)
    paramat <- do.call('rbind', params)
    paramat <- data.frame(apply(paramat, 2, as.numeric))
    if (distr == 'zinegbin') {
      data <- rmvzinegbin(n, ks=paramat$size, munbs=paramat$munb,
                  ps=paramat$pstr0, Sigma=Sigma)
    }
    structure(list(comm=data, params=if (retParams) params else NULL), class = "community")
}

#' Get the parameters for the OTUs (rows) of each community
#'
#' @param comm community: matrix of counts
#' @param distr distribution to fit (see fitdistr)
#' @return list of parameters
#' @export
get_comm_params <- function(comm, distr, ...) {

    lapply(1:nrow(comm),  function(i) {
        x <- as.numeric(comm[i,])
        ll <- c(list(), fitdistr(x, distr, ...)$par)
        ll$mean <- mean(x)
        ll
    })
}

#' Fit parameters of a marginal distribution to some data vector
#'
#' @param x data vector
#' @param densfun string giving distribution function name
#' @param start starting guess for the parameters (recommended leaving this out)
#' @param ... further arguments to densfun
#' @importFrom VGAM dzinegbin dzipois
#' @importFrom stats dnbinom
#' @export
fitdistr <- function (x, densfun, start, ...)  {
    if (class(x) != "numeric") stop("Error: input must be numeric vector")
    Call <- match.call(expand.dots = TRUE)
    if (missing(start)) 
        start <- NULL
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]
    if (missing(x) || length(x) == 0L || mode(x) != "numeric") 
        stop("'x' must be a non-empty numeric vector")
    if (any(!is.finite(x))) 
        stop("'x' contains missing or infinite values")
    if (missing(densfun) || !(is.function(densfun) || is.character(densfun))) 
        stop("'densfun' must be supplied as a function or name")
    control <- list()
    n <- length(x)
    if (is.character(densfun)) {
        distname <- tolower(densfun)
        densfun <- switch(distname, zipois = VGAM::dzipois, zinegbin = VGAM::dzinegbin, negbin = dnbinom,
            pois=dpois, NULL)
        if (is.null(densfun)) 
            stop("unsupported distribution")
    }
    if (distname == "zipois") {
        if (!is.null(start)) 
            stop(gettextf("supplying pars for the %s distribution is not supported", 
              "Poisson"), domain = NA)
        whichz  <- which(x == 0.0)
        which1  <- which(x == 1.0)
        max <- abs(length(whichz) - length(which1))
        max <- max - max*.1
        zind    <- whichz[1:max]
        tempx   <- x[-zind]
        pstr0  <- length(which(x == 0)) / length(x)
        pstr0  <- ifelse(pstr0 == 0, 0, abs(pstr0 - exp(-mean(tempx))))   # correct for approx expected zeros in a poisson (important for small rates)
        estimate <- mean(x) / (1 - pstr0)
        vars <- ((1 - pstr0) * (estimate^2 + estimate)) - ((1-pstr0) * estimate)^2
        sds  <- sqrt(vars)
        start <- c(lambda = estimate, pstr0 = pstr0)
        start <- start[!is.element(names(start), dots)]
        lower <- c(1e-4, 0)
        upper <- c(Inf, .99)
        loglikfn <- match.fun(logLikzip)
        }
    if (distname == "gamma" && is.null(start)) {
        if (any(x < 0)) 
            stop("gamma values must be >= 0")
        m <- mean(x)
        v <- var(x)
        start <- list(shape = m^2/v, rate = m/v)
        start <- start[!is.element(names(start), dots)]
        control <- list(parscale = c(1, start$rate))
    }
    else if (distname == "zinegbin" && is.null(start)) {
        whichz  <- which(x == 0.0)
        which1  <- which(x == 1.0)
        max   <- abs(length(whichz) - length(which1))
        max   <- max - max*.1
        zind  <- na.omit(whichz[1:max])
        tempx <- x[-zind]
        pstr0 <- length(which(x == 0)) / length(x)
        pstr0 <- ifelse(pstr0 == 0, 0, abs(pstr0 - exp(-mean(tempx))))   # correct for approx expected zeros, if there are any
        m     <- mean(x)
        estimate <- m / (1 - pstr0)
        v    <- var(x)
        if (v < m) size <- 1e-3
        else size <- estimate / ((v/m) - 1 - (pstr0 * estimate))
        start <- c(size = size, munb = estimate, pstr0 = pstr0)
        start <- start[!is.element(names(start), dots)]
        lower <- c(1e-4, 1e-4, 0)
        upper <- c(1e4, 1e4, .99)
        loglikfn <- match.fun(logLikzinb)
     } 
    else if (distname == 'negbin') {
        m <- mean(x)
        v <- var(x)
        size <- if (v > m) 
                m^2/(v - m)
                else 100
        start <- list(size = size, mu = m)
        start <- start[!is.element(names(start), dots)]
        lower <- c(1e-4, 1e-4)
        upper <- c(1e4, 1e4)
        loglikfn <- match.fun(logLiknb)
    }
    else if (distname == "pois") {
        m <- mean(x)
    return(list(par=list(lambda=m)))
    }

    res <- optim(start, loglikfn, x=x, method='L-BFGS-B', lower=lower, upper=upper, ...)
    return(res)
}


#' @keywords internal
#' @importFrom VGAM dzinegbin
logLikzinb <- function(param, x, ddistr, ...) {
    pstr0  <- param['pstr0']
    munb <- param['munb']
    size   <- param['size']
    -sum(VGAM::dzinegbin(x, munb=munb, pstr0=pstr0, size=size, log=TRUE, ...))
}

#' @keywords internal
#' @importFrom stats dnbinom
logLiknb <- function(param, x, ddistr, ...) {
    munb <- param['mu']
    size   <- param['size']
    -sum(dnbinom(x, mu=munb, size=size, log=TRUE, ...))
}
#' @keywords internal
#' @importFrom VGAM dzipois
logLikzip <- function(param, x, ddistr, ...) {
    pstr0  <- param['pstr0']
    lambda <- param['lambda']
    -sum(VGAM::dzipois(x, lambda=lambda, pstr0=pstr0, log=TRUE, ...))
}

#' qq-plot for theoretical vs observed communities
#'
#' @export
qqdplot_comm <- function(comm, distr, param, plot=TRUE, ...) {
    if (class(comm) != 'qqdcomm') {
        if (!missing(param)) {
            fit <- t(sapply(1:nrow(comm), function(i) qqdplot(comm[i,], distr=distr, param=param[[i]], plot=FALSE)))
        } else {
            fit <- t(apply(comm, 1, qqdplot, distr=distr, param=NULL, plot=FALSE))
        }
        x   <- as.vector(fit)
        y   <- as.vector(as.matrix(comm))
    } else {
        x <- comm$x
        y <- comm$y
    }
    
    if (plot) {
        xrange <- range(y)
        if (!('main' %in% names(list(...)))) main = 'QQ-plot'
        else main <- list(...)$main
        plot(y, x, main=main, ylim=xrange, xlab='Observed Quantiles', ylab='Theoretical Quantiles')
        abline(0, 1, lty=2)
        legend('topleft', paste('r^2 =', format(summary(lm(x ~ y))$adj.r.squared, digits=4), sep=""), box.lwd=0)
    } else {
        return(structure(list(x=x, y=y), class="qqdcomm"))
    }
}

#' deprecated??
#' @keywords internal
qqdplot <- function(y, distr, param, plot=TRUE, ...) {

    if ((missing(param) || is.null(param)) && !(distr %in% c('lnorm', 'pois', 'negbin'))) {
        param <- as.list(fitdistr(y, distr, ...)$par)
    }
    y <- as.numeric(y)
    n <- length(y)
    if (distr == 'zinegbin') {
        x <- qzinegbin(ppoints(n), munb=param$munb, size=param$size, pstr0=param$pstr0)[order(order(y))]
    } else if (distr == 'zipois' ) {
        x <- qzipois(ppoints(n), lambda=param$lambda, pstr0=param$pstr0)[order(order(y))]
    } else if (distr == 'lnorm') {
        logy <- log1p(y)
        sd   <- sd(logy)
        mean <- mean(logy)
        x <- qlnorm(ppoints(n), meanlog=mean, sdlog=sd)[order(order(y))]
    } else if (distr == 'negbin') {
        if (missing(param) || is.null(param)) {
            require(MASS)
            param <- as.list(fitdistr(y, 'negbin')$par)
        }
        x <- qnbinom(ppoints(n), mu=param$mu, size=param$size)[order(order(y))]
    } else if (distr == 'pois' || is.null(param)) {
        if (missing(param) || is.null(param)) {
            require(MASS)
            param <- as.list(fitdistr(y, 'poisson')$estimate)
        }
        x <- qpois(ppoints(n), lambda=param$lambda)[order(order(y))]
    }
    
    if (plot) {
        xrange <- range(y)
        plot(y, x, main="QQ-plot", ylim=xrange, xlab='Observed Quantiles', ylab='Theoretical Quantiles')
        abline(0, 1, lty=2)
        legend('topleft', paste('r^2 =', format(summary(lm(x ~ y))$adj.r.squared, digits=4), sep=""), box.lwd=0)
    } else {
        return(x)
    }
}



