#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ReflectIndexMed(NumericVector vec, int ind, std::string side) {
  int n = vec.size();
  double out = NA_REAL;
  double med;
  int dist_to_end;
  int dist_to_go;
  double median_dist;
  if (ind >= 0 && ind < n) {
    if (side == "left") {
      NumericVector left(ind);
      for (int i = 0; i < ind; i++) {
        left[i] = vec[i];
      }
      med = median(left);
      dist_to_end = ind;
      dist_to_go = 2 * dist_to_end;
      median_dist = 1 + 0.5 * (dist_to_end - 1);
      out = vec[ind] + dist_to_go / median_dist * (med - vec[ind]);
    } else if (side == "right") {
      NumericVector right(n - (ind + 1));
      for (int i = 0; i < n - (ind + 1); i++) {
        right[i] = vec[ind + i + 1];
      }
      med = median(right);
      dist_to_end = (n - 1) - ind;
      dist_to_go = 2 * dist_to_end;
      median_dist = 1 + 0.5 * (dist_to_end - 1);
      out = vec[ind] + dist_to_go / median_dist * (med - vec[ind]);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector Smooth(NumericVector vec) {
  int n = vec.size();
  NumericVector smoothed(n);
  if (n > 1) {
    NumericVector three = vec[IntegerVector::create(0, 0, 1)];
    smoothed[0] = mean(three);
    three = vec[IntegerVector::create(n - 2, n - 1, n - 1)];
    smoothed[n - 1] = mean(three);
    for (int i = 1; i < n - 1; i++) {
      three = vec[IntegerVector::create(i - 1, i, i + 1)];
      smoothed[i] = mean(three);
    }
  } else {
    smoothed[0] = vec[0];
  }
  return smoothed;
}

// [[Rcpp::export]]
NumericVector MedReflectExtend(NumericVector vec, bool preserve_mean = false,
                               bool smooth = false) {
  int n = vec.size();
  if (n <= 1) {
    return vec;
  } else {
    NumericVector extended((n - 1) + n + (n - 1));
    for (int i = 0; i < n - 1; i++) {
      extended[(n - 2) - i] = ReflectIndexMed(vec, i + 1, "left");
    }
    for (int i = 0; i < n; i++) {
      extended[i + (n - 1)] = vec[i];
    }
    for (int i = 0; i < n - 1; i++) {
      extended[(n - 1) + n + i] = ReflectIndexMed(vec, (n - 1) - 1 - i,
                                                  "right");
    }
    if (smooth) {
      NumericVector side = extended[seq(0, n - 2)];
      NumericVector smoothed = Smooth(side);
      for (int i = 0; i < n - 1; i++) {
        extended[i] = smoothed[i];
      }
      side = extended[seq((n - 1) + n, (n - 1) + n + (n - 2))];
      smoothed = Smooth(side);
      for (int i = 0; i < n - 1; i++) {
        extended[(n - 1) + n + i] = smoothed[i];
      }
    }
    if (preserve_mean) {
      double mean_vec = mean(vec);
      double to_add;
      int half_n_to_add_to = n - 1 - 1;
      int atoms_to_add = 2 * half_n_to_add_to * (half_n_to_add_to + 1) / 2;
      double atom;
      int billion = pow(10, 9);
      while (fabs(mean_vec - mean(extended)) > mean_vec / billion) {
        // this procedure suffers from some numerical imprecision
        // but repeating it many times yields the required result
        to_add = mean(vec) * extended.size() - sum(extended);
        atom = to_add / atoms_to_add;
        for (int i = 1; i < n - 1; i++) {
          extended[(n - 2) - i] += i * atom;
        }
        for (int i = 1; i < n - 1; i++) {
          extended[(n - 1) + n + i] += i * atom;
        }
      }
    }
    return extended;
  }
}

// [[Rcpp::export]]
NumericMatrix MedReflectExtendRows(NumericMatrix rows,
                                   bool preserve_mean = false,
                                   bool smooth = false) {
  NumericMatrix extended_rows(rows.nrow(), 3 * rows.ncol() - 2);
  for (int i = 0; i < rows.nrow(); i++) {
    extended_rows(i, _) = MedReflectExtend(rows(i, _), preserve_mean, smooth);
  }
  return extended_rows;
}

//' Exponentially smooth a series of observations.
//'
//' This function assumes that the observations are evenly spaced and separated
//' by 1 time unit (so choose your \code{tau} based on that).
//'
//' @param obs A numeric vector of observations (in order).
//' @param tau The time scale for the exponential smoothing (see Stroud 1999).
//' @param extended Logical. Has the series (`obs`) already been extended via
//'   via `nandb:::MedReflectExtend()`? If not, \code{ExpSmooth} will do this
//'   prior to smoothing as an edge-correction technique. You will probably
//'   never set this to \code{TRUE}, but [BestTau()] needs this feature.
//'
//' @return The smoothed series, a numeric vector of the same length.
//' @examples
//' ExpSmooth(1:10, 1)
//' @export
// [[Rcpp::export]]
NumericVector ExpSmooth(NumericVector obs, double tau, bool extended = false) {
  int n = obs.size();
  double numerator;
  double denominator;
  NumericVector obs_extended;
  if (extended) {
    if ((n + 2) % 3 != 0) {
      return rep(NA_REAL, n);
    }
    obs_extended = obs;
    n = (n + 2) / 3;
  } else {
    obs_extended = MedReflectExtend(obs, true, true);
  }
  NumericVector smoothed(n);
  NumericVector weights(2 * n - 1);
  for (int i = 0; i < 2 * n - 1; i++) {
    weights[i] = exp(- i / tau);
  }
  int m = obs_extended.size();
  int i_in_extended;
  for (int i = 0; i < n; i++) {
    numerator = 0;
    denominator = 0;
    for (int j = 0; j < m; j++) {
      i_in_extended = (n - 1) + i;
      numerator += weights[abs(j - i_in_extended)] * obs_extended[j];
      denominator += weights[abs(j - i_in_extended)];
    }
    smoothed[i] = numerator / denominator;
  }
  return smoothed;
}

// Exponentially smooth a series of observations without edge correction.
//
// This function assumes that the observations are evenly spaced and separated
// by 1 time unit (so choose your \code{tau} based on that).
//
// The function `ExpSmooth()` performs edge correction. This one does not.
//
// @param obs A numeric vector of observations (in order).
// @param tau The time scale for the exponential smoothing (see Stroud 1999).
//
// @return The smoothed series, a numeric vector of the same length.
// @examples
// ExpSmoothNaive(1:10, 1)
// [[Rcpp::export]]
NumericVector ExpSmoothNaive(NumericVector obs, double tau) {
  int n = obs.size();
  double numerator;
  double denominator;
  NumericVector smoothed(n);
  NumericVector weights(n);
  for (int i = 0; i < n; i++) {
    weights[i] = exp(- i / tau);
  }
  for (int i = 0; i < n; i++) {
    numerator = 0;
    denominator = 0;
    for (int j = 0; j < n; j++) {
      numerator += weights[abs(j - i)] * obs[j];
      denominator += weights[abs(j - i)];
    }
    smoothed[i] = numerator / denominator;
  }
  return smoothed;
}

//' Exponentially smooth pillars of a 3-dimensional array
//'
//' For a 3-dimensional array \code{mat3d}, pillar \code{i,j} is defined as
//' \code{mat3d[i, j, ]}. \code{ExpSmoothPillars} function performs
//' \link{ExpSmooth} on each pillar. \code{ExpSmoothRows} performs
//' \link{ExpSmooth} on each row of a matrix.
//'
//' @param mat3d A 3-dimensional array.
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
  if (R_IsNA(tau)) {
    return smoothed_pillars;
  }
  IntegerVector dim = mat3d.attr("dim");
  int n_pillars = dim[0] * dim[1];
  int pillar_len = dim[2];
  NumericVector pillar_i(pillar_len);
  NumericVector smoothed_pillar_i(pillar_len);
  for (int i = 0; i < n_pillars; i++) {
    for (int j = 0; j < pillar_len; j++) {
      pillar_i[j] = mat3d[i + j * n_pillars];
    }
    smoothed_pillar_i = ExpSmooth(pillar_i, tau);
    for (int j = 0; j < pillar_len; j++) {
      smoothed_pillars[i + j * n_pillars] = smoothed_pillar_i[j];
    }
  }
  return smoothed_pillars;
}

//' @rdname ExpSmoothPillars
//' @param mat A matrix.
//' @param extended Logical. Has the series (`obs`) already been extended via
//'   via `nandb:::MedReflectExtend()`? If not, \code{ExpSmooth} will do this
//'   prior to smoothing as an edge-correction technique. You will probably
//'   never set this to \code{TRUE}, but [BestTau()] needs this feature.
//' @export
// [[Rcpp::export]]
NumericMatrix ExpSmoothRows(NumericMatrix mat, double tau,
                            bool extended = false) {
  if (extended) {
    NumericMatrix smoothed_rows(mat.nrow(), (mat.ncol() + 2) / 3);
    for (int i = 0; i < mat.nrow(); i++) {
      smoothed_rows(i, _) = ExpSmooth(mat(i, _), tau, true);
    }
    return smoothed_rows;
  } else {
    NumericMatrix smoothed_rows(clone(mat));
    for (int i = 0; i < mat.nrow(); i++) {
      smoothed_rows(i, _) = ExpSmooth(mat(i, _), tau);
    }
    return smoothed_rows;
  }
}
