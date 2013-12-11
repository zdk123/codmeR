#!/usr/bin/env Rscript

##########################
# @author Zachary Kurtz
#
#  Simulate communities from log-normal, poisson or negative binomial distributions
######

log_norm_mean_vector <- function(D, logmean=1, logsd=.5) {
# Generate an micrbial community-like vector of D means
    rlnorm(D, logmean,logsd)
}

fixed_mean_vector <- function(D, neff, val=2, numhigher=ceiling(D/50)) {
    if (neff > D) stop("Effective species number is higher than number of species")
    y <- optim(rep(val, numhigher), fitfn, x=(x<-rep(val, D-numhigher)), neff=neff, method="L-BFGS-B", lower=val, upper=Inf)
    sample(c(x,y$par))
}

fitfn <- function(y, x, neff) {
    x.f <- c(x,y)/sum(c(x,y))
    x <- x.f[1:length(x)]
    y <- x.f[(length(x)+1):length(x.f)]
    abs(log(neff) + sum(x*log(x)) + sum(y*log(y)))*1e4
}

sim_counts <- function(n = 1, mu, Sigma=NULL, method='ZIP', negscheme=exp, targetneff, ...) {
# Simulate count data 
# Args:
#   n -> number of samples
#   mu -> vector of 'OTU' means across samples
#   Sigma -> positive definite matrix, one will be generated if not provided
#   negscheme -> function to deal with negative values
#   required: length(mu) == nrow(Sigma) == ncol(Sigma)

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


unsparseSigma <- function(dimension, avecor=0.5, sparsity=0.5) {

    require(msm) ; require(Matrix)
    cormat <- matrix(0, dimension, dimension)
    d <- dimension

    cors <- rtnorm((d*(d-1))/2, mean=avecor, lower=-.999, upper=.999, sd=sparsity)
    cormat[lower.tri(cormat, diag=F)] <- cormat[upper.tri(cormat, diag=F)] <- cors
    diag(cormat) <- 1
    return(nearPD(cormat, corr=TRUE)$mat)
}


addmin <- function(x, jitter=1e-3) {

    x + abs(min(x)) + jitter

}


shannon <- function(x) {
    x.f <- (x+1)/sum(x+1)
    -sum(x.f*log(x.f))
}

neff <- function(x) exp(shannon(x))

sparseSigma <- function(dimension, sparsity=0.6, center=0.0) {
#Generates a sparse, positive definite [square] correlation matrix
#Args:
#   dimension -> number of rows/cols
#   sparsity -> controls probability of sparsity by setting the standard deviation 
#                of truncated normal distribution. Lower SD -> higher probability of a draw
#                closer to the center.

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

sparseSigma_JF <- function(D, p, var=1, center=0.0) {
# Jonathan Friedman's described method of computing a correlation matrx
# Args:
# D -> Dimension
# p -> probability that any otu pairs are perfectly correlated
# var -> individual variance (same for all OTUs)
# center - > average correlation (default 0.0)
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


convert_MD_format <- function(matrix, textname="otus") {
    if (is.null(colnames(matrix))) colnames(matrix) <- 1:ncol(matrix)
    if (is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
    header <- c(textname, colnames(matrix))
    newmat <- cbind(rownames(matrix), matrix)
    newmat <- rbind(header, newmat)
    return(newmat)
}


main <- function() {
suppressPackageStartupMessages(require(optparse))

#define option list, -h by default
option_list <- list(
    make_option(c("-d", "--num_OTUs"), action="store", default=8, type="numeric",
        help="Number of OTUs [default: %default]"),
    make_option(c("-n", "--num_Samples"), action="store", type="integer", default=30,
        help="Number of Samples [default: %default]"),
    make_option(c("-o", "--data_out_fp"), action="store", type="character",
        help="file path to write data [no default]"),
    make_option(c("-c", "--cor_out_fp"), action="store", type="character",
        help="file path & basename to write input correlation matrix [no default]"),
    make_option(c("-a", "--ave_corr"), action="store", default=0,
        help="Average value of correlations [default %default]"),
    make_option(c("-s", "--sparse"), action="store_true", default=FALSE,
        help="Use SMVLNparse method to generate correlation matrix (most correlations are zero, spread is controlled by sd_param, center is ave_corr [default: %default]) "),
    make_option(c("-p", "--sparse_param"), action="store", default=2, type='numeric',
        help="Sparsity (sd of correlations) or s*D target edges (if -qG) [default: %default]"),
    make_option(c("-i", "--cov_stats_fp"), action="store", default="",
        help="file path to store covariance summary statistics [default: Don't Save]"),
    make_option(c("-t", "--mean_method"), action="store", default="log-norm",
        help="Pick method for drawing mean vector, log-norm (with sd parameter) or fixed (target neff) [default %default]"),
    make_option(c("-v", "--div_param"), action="store", default=0.5, type="numeric",
        help="Standard deviation of mean if log norm (0-.8) Neff if fixed [default: %default]"),
    make_option(c("-G", "--graph_method"), action="store", default="",
        help="Use graph methods to create graph (new/improved & prefered)"),
    make_option(c("-m", "--method"), action="store", default="ZIP",
        help="method for simulating data. MVLN for multivariate log-normal, MN - multinomial, ZIP - zero-inflated poisson, ZINB - zero-inflated negative binomial [default: %default]"),
    make_option(c("-N", "--normalize"), action="store_true", default=FALSE,
        help="Normalize results? Perfoms a variety of procedures [Default: \"\"]"),
    make_option(c("-l", "--parallel"), action="store_true", default=FALSE,
        help="parallelize methods where possible")
    )
    opt <- parse_args(OptionParser(option_list=option_list))
  # run with defaults
    args<-commandArgs(TRUE)
    D <- as.numeric(opt$num_OTUs)
    n <- as.numeric(opt$num_Samples)
    sd <- as.numeric(opt$sparse_param)
    if (opt$graph_method != "") {
        if (!(opt$graph_method %in% c("erdos_renyi", "cluster", "band", "scale_free", "hub"))) {
             stop("Graph method not currently supported")
        }
        source("./community_graph.R")
        require(Matrix)
        #target sparsity pattern O(1/n)
        constS <- opt$sparse_param
        sp     <- constS*D/((D-1)*D/2)
        edges  <- sp*(((D-1)*D/2))
        
        Cov <- matrix(-1, D, D)
        it  <- 1
        while (any(sign(diag(Cov)) == -1)) {
            Graph <- make_graph(opt$graph_method, D, e=edges)
            Prec  <- graph_to_precision(Graph)
            Cov  <- prec2cov(Prec)
            if (it > 50) {
                stop("Error: covariance structure is not positive definite")
                break
            }
            it <- it + 1
        }
        
        Sigma <- cov2cor(Cov)
        
        
    } else {
        if (opt$sparse)     Sigma <- sparseSigma_JF(D, sd, center = as.numeric(opt$ave_corr))
        else Sigma <- unsparseSigma(D, as.numeric(opt$ave_corr), sd)
    }

    if (opt$method %in% c('ZIP', 'ZINB')) meanOfMu <- 4
    else meanOfMu <- 0.3
    
    if (opt$mean_method == "log-norm") {
        mu <- log_norm_mean_vector(D, meanOfMu, opt$div_param)
    } else if (opt$mean_method == "fixed") {
        if (opt$div_param <= 1 || opt$div_param > D) stop("Error: invalid target neff")
        mu <- fixed_mean_vector(D, val=meanOfMu, neff=opt$div_param)
    } else {
        stop("mean vector method is not supported")
    }

    data <- sim_counts(n, mu, Sigma, method=opt$method, parallel=opt$parallel)
    Y <- t(data$data)
    
    cor <- data$empcor
    trucor <- data$incor
    fout <- strsplit(opt$cor_out_fp, "\\.")[[1]]
    if (length(fout) == 1) fout[2] <- "txt"
    empcorfile <- paste(fout[1], "_emp.", fout[2], sep="")
    trucorfile <- paste(fout[1], "_input.", fout[2], sep="")


    rownames(Y) <- 1:nrow(Y)
    colnames(Y) <- 1:ncol(Y)
    rownames(cor) <- 1:nrow(cor)
    colnames(cor) <- 1:ncol(cor)
    rownames(trucor) <- 1:nrow(trucor)
    colnames(trucor)   <- 1:ncol(trucor)

    newY <- convert_MD_format(Y)
    newcor <- convert_MD_format(cor)
    newtru <- convert_MD_format(trucor)

    write.table(newY, opt$data_out_fp, sep="\t", col.names=F, row.names=F, quote=F)
#    write.table(newcor, empcorfile, sep="\t", col.names=F, row.names=F, quote=F)
    write.table(newtru, trucorfile, sep="\t", col.names=F, row.names=F, quote=F)
    
    if (opt$normalize != "") {
        source("./normalization.R")
        Ynormlist <- list(
        Y.log  = apply(Y, 2, logpseudo),
        Y.tsum = (Y.tsum <- apply(Y, 2, norm_to_total)),
        Y.psum = (Y.psum <- apply(Y, 2, norm_pseudo)),
        Y.diri = (Y.diri <- apply(Y, 2, norm_diric)),
        Y.clr  = clr(Y.tsum),
        Y.clrp = clr(Y.psum),
        Y.clrd = clr(Y.diri),
        Y.CSS  = CSS(Y),
        Y.DES  = DESeq(Y)
        )
          
       dfout <- strsplit(opt$data_out_fp, "\\.")[[1]]
       if (length(dfout) == 1) dfout[2] <- "txt"
       dfoutlist <- list("_LogCounts.", "_TotalSumNorm.", "_TotalSumPseudo.", 
                          "_DirichletNorm.", "_CLR.", "_PseudoCLR.", 
                          "_DirichletCLR.", "_CSS.", "_DESeq.")
       lapply(1:length(Ynormlist), function(i) {
            write.table(convert_MD_format(Ynormlist[[i]]), paste(dfout[1], dfoutlist[[i]], dfout[2], sep=""), 
                sep="\t", col.names=F, row.names=F, quote=F)
       })
    }
    
    
    if (opt$cov_stats_fp != "") {
        source("./community_graph.R")
        gR <- graphReport(Graph)
        cR <- covReport(Cov, Prec)
        reports <- c(gR, cR)
        sink <- lapply(1:length(reports), function(i) {
                    name <- names(reports)[i]
                    tempreport <- as.matrix(reports[[i]])
                    colnames(tempreport) <- rep("", ncol(tempreport))
                    if (i != length(reports)) rownames(tempreport) <- rep("", nrow(tempreport))
                    data <- convert_MD_format(tempreport, textname=name)
                    write(data, file=opt$cov_stats_fp, append=TRUE, ncolumns=2000) })
    }
    
}



if (!interactive()) {
#options(warn=-1)
suppressPackageStartupMessages(main())
}



cor2beta <- function(cor, data) {
    sds <- apply(data, 2, sd)
    cor * (sds %*% t(1/sds))
}


