#include <Rcpp.h>
using namespace Rcpp;

//' Exponentially smooth a series of observations.
//'
//' This function assumes that the observations are evenly spaced and separated
//' by 1 time unit (so choose your \code{tau} based on that).
//'
//' @param obs A numeric vector of observations (in order).
//' @param tau The time scale for the exponential smoothing (see Stroud 1999).
//'
//' @return The smoothed series, a numeric vector of the same length.
//' @examples
//' ExpSmooth(1:10, 1)
//' @export
// [[Rcpp::export]]
NumericVector ExpSmooth(NumericVector obs, double tau) {
  int n = obs.size();
  NumericVector weights0(n);
  for (int i = 0; i < n; i++) {
    weights0[i] = exp(- i / tau);
  }
  NumericVector weights(n);
  NumericVector new_obs(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      weights[j] = weights0[abs(i - j)];
    }
    new_obs[i] = sum(weights * obs) / sum(weights);
  }
  return new_obs;
}

//' Exponentially smooth pillars of a 3-dimensional array
//'
//' For a 3-dimensional array \code{mat3d}, pillar \code{i,j} is defined as
//' \code{mat3d[i, j, ]}. \code{ExpSmoothPillars} function performs
//' \link{ExpSmooth} on each pillar. \code{ExpSmoothRows} performs
//' \link{ExpSmooth} on each row of a matrix.
//'
//' @param mat3d A 3-dimensional array.
//' @param mat A matrix.
//' @param tau The time scale for the exponential smoothing (see Stroud 1999).
//'
//' @return For \code{ExpSmoothPillars}, a 3-dimensional array where each
//' pillar has been smoothed. For \code{ExpSmoothRows}, a matrix where each
//' row has been smoothed.
//'
//' @examples
//' m3 <- array(1:12, dim = c(2, 2, 3))
//' ExpSmoothPillars(m3, 7)
//' @export
// [[Rcpp::export]]
NumericVector ExpSmoothPillars(NumericVector mat3d, double tau) {
  // mat3d is actually passed in as a 1d vector
  NumericVector smoothed_pillars(clone(mat3d));
  IntegerVector dim = mat3d.attr("dim");
  int n_pillars = dim[0] * dim[1];
  int pillar_len = dim[2];
  NumericVector pillar_i(pillar_len);
  for (int i = 0; i < n_pillars; i++) {
    for (int j = 0; j < pillar_len; j++) {
      pillar_i[j] = mat3d[i + j * n_pillars];
    }
    pillar_i = ExpSmooth(pillar_i, tau);
    for (int j = 0; j < pillar_len; j++) {
      smoothed_pillars[i + j * n_pillars] = pillar_i[j];
    }
  }
  return smoothed_pillars;
}

//' @rdname ExpSmoothPillars
//' @export
// [[Rcpp::export]]
NumericMatrix ExpSmoothRows(NumericMatrix mat, double tau) {
  NumericMatrix smoothed_rows(clone(mat));
  for (int i = 0; i < mat.nrow(); i++) {
    smoothed_rows(i, _) = ExpSmooth(mat(i, _), tau);
  }
  return smoothed_rows;
}
