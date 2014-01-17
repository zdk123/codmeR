#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
 
arma::mat distsq(arma::mat A) {
    int m = A.n_rows, 
        k = A.n_cols;
    arma::colvec An =  sum(square(A),1);
    arma::mat C = -2 * (A * A.t());
    C.each_col() += An;
    C.each_row() += An.t();
    return C; 
}

arma::mat sumdiff(arma::mat A) {
    int m = A.n_rows, 
        k = A.n_cols;
    arma::mat onemat = arma::ones(k, m);
    arma::mat rowsums = A * onemat;
    return rowsums.t() - rowsums;
}

// [[Rcpp::export]]
arma::mat fastaitvar(arma::mat A) {
    int m = A.n_rows, 
        k = A.n_cols;
    arma::mat Alog = log(A);
    return (distsq(Alog) - square(sumdiff(Alog))/k) / (k-1);
}
