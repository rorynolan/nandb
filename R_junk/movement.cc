#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int MostConsecutiveLEs(NumericVector x, double thresh) {
  if (any(is_na(x)))
    return NA_INTEGER;
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

// [[Rcpp::export]]
IntegerMatrix MostConsecutiveLEsPillars(NumericVector mat3d, double thresh) {
  IntegerVector dim = mat3d.attr("dim");
  int n_pillars = dim[0] * dim[1], pillar_len = dim[2];
  IntegerMatrix mcles(dim[0], dim[1]);
  NumericVector pillar_i(pillar_len);
  for (int i = 0; i < n_pillars; i++) {
    for (int j = 0; j < pillar_len; j++) {
      pillar_i[j] = mat3d[i + j * n_pillars];
    }
    mcles(i % dim[0], i / dim[0]) = MostConsecutiveLEs(pillar_i, thresh);
  }
  return mcles;
}

