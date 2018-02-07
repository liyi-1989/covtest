#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List trand(NumericVector x) {
  
  int n=100;
  NumericVector v=rnorm(n);
  double M=mean(v);
  double S=sd(v);
  return List::create(M,S);
}

