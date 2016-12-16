#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector WhichIntervalC(NumericVector numbers, NumericMatrix ranges) {
  int nn = numbers.size();
  IntegerVector interval(nn);
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < ranges.nrow(); j++) {
      if (numbers[i] > ranges(j, 0) && numbers[i] <= ranges(j, 1)) {
        interval[i] = j + 1;
        break;
      }
    }
  }
  interval.attr("dim") = numbers.attr("dim");
  return(interval);
}
