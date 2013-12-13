#include <Rcpp.h>
using namespace Rcpp;

//' compute aitchison variation of matrix cols quickly
//' @param x NumericMatrix with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
NumericMatrix fastaitvar(NumericMatrix x) {
    int m = x.nrow(), 
        k = x.ncol();
    NumericMatrix C(k, k);

    for (int i=0; i!=m; ++i)
      {
      for (int j=(i+1); j!=m; ++j)
        {
            C(i,j) = var(log(x.column(i) / x.column(j)));
            C(j,i) = C(i,j);
        }
      }
    return C;
}

