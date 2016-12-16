#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int MostConsecutiveLEs(NumericVector x, double thresh) {
  LogicalVector y = x <= thresh;
  int i = 0, j = 0, most = 0;
  while (i < x.size()) {
    if (y[i])
      j++;
    else {
      if (j > most)
        most = j;
      j = 0;
    }
    i++;
  }
  if (j > most)
    most = j;
  return most;
}


