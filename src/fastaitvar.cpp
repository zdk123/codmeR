#include <Rcpp.h>
using namespace Rcpp;

//' compute aitchison variation vector
//' @param x NumericMatrix with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
double T(NumericVector xi, NumericVector xj) {
    return var(log(xi / xj));
}

//' compute aitchison variation of matrix cols quickly
//' @param x NumericMatrix with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
NumericMatrix fastaitvar(NumericMatrix x) {
    int m = x.nrow(), 
        k = x.ncol();
    Rcpp::NumericMatrix C(k, k);

    for (int i=0; i!=k; ++i)
      {
      for (int j=(i+1); j!=k; ++j)
        {
            C(i,j) = T(x.column(i), x.column(j));
            C(j,i) = C(i,j);
        }
      }
    return C;
}

