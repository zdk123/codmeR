#include <Rcpp.h>
using namespace Rcpp;

//' compute Hamming distance
//' @param xi NumericVector with data
//' @param xj NumericVector with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
int Hamdist(NumericVector xi, NumericVector xj) {
    return sum(xi != xj);
}

//' compute Hamming distance of matrix cols quickly
//' @param x NumericMatrix with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
NumericMatrix HamdistMat(NumericMatrix x) {
    int m = x.nrow(),
        k = x.ncol();
    Rcpp::NumericMatrix C(k, k);

    for (int i=0; i!=k; ++i)
      {
      for (int j=(i+1); j!=k; ++j)
        {
            C(i,j) = Hamdist(x.column(i), x.column(j));
            C(j,i) = C(i,j);
        }
      }
    return C;
}


// [[Rcpp::export]]
int coCalc(NumericVector xi, NumericVector xj, int sumint = 2) {
    return sum((sign(xi) + sign(xj)) == sumint);
}

// [[Rcpp::export]]
NumericMatrix coCalcMat(NumericMatrix x, int sumint = 2) {
    int m = x.nrow(),
        k = x.ncol();
    Rcpp::NumericMatrix C(k, k);

    for (int i=0; i!=k; ++i)
      {
      for (int j=i; j!=k; ++j)
        {
            C(i,j) = coCalc(x.column(i), x.column(j), sumint);
            C(j,i) = C(i,j);
        }
      }
    return C;
}

//' compute Cooccurance
//' @param xi NumericVector with data
//' @param xj NumericVector with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
NumericMatrix cooccurMat(NumericMatrix x) {
    return coCalcMat(x, 2);
}

//' compute Coabsence
//' @param xi NumericVector with data
//' @param xj NumericVector with data
//' @return symmetric matrix with pairwise Aitchison variations
// [[Rcpp::export]]
NumericMatrix coabsenceMat(NumericMatrix x) {
    return coCalcMat(x, 0);
}



