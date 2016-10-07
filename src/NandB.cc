// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/sign.hpp>

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
//' \code{mat3d[i, j, ]}. This function performs \link{ExpSmooth} on each pillar.
//'
//'
//' @param tau The time scale for the exponential smoothing (see Stroud 1999).
//'
//' @return A 3-dimensional array where each pillar has been smoothed.
//'
//' @examples
//' m3 <- array(1:12, dim = c(2, 2, 3))
//' ExpSmoothPillars(m3, 7)
//' @export
// [[Rcpp::export]]
NumericVector ExpSmoothPillars(NumericVector mat3d, double tau) {
	// mat3d is actually passed in as a 1d vector
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
			mat3d[i + j * n_pillars] = pillar_i[j];
		}
	}
	mat3d.attr("dim") = dim;
	return mat3d;
}

//' Get the means/medians/variances of pillars of a 3d array
//'
//' For a 3-dimensional array \code{mat3d}, pillar \code{ij} is defined as
//' \code{mat3d[i, j, ]}. These functions compute the mean, median and variance of each
//' pillar.
//'
//' @param mat3d A 3-dimensional array.
//'
//' @return A matrix where element \code{i,j} is equal to \code{mean(mat3d[i, j, ])},
//' \code{median(mat3d[i, j, ])}, or \code{var(mat3d[i, j, ])}.
//'
//' @examples
//' m3 <- array(1:16, dim = c(2, 2, 4))
//' MeanPillars(m3)
//' MedianPillars(m3)
//' VarPillars(m3)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix MeanPillars(NumericVector mat3d) {
  IntegerVector dim = mat3d.attr("dim");
	int n_pillars = dim[0] * dim[1];
	int pillar_len = dim[2];
	NumericMatrix means(dim[0], dim[1]);
	double mean_i;
	NumericVector pillar_i(pillar_len);
	for (int i = 0; i < n_pillars; i++) {
		for (int j = 0; j < pillar_len; j++) {
			pillar_i[j] = mat3d[i + j * n_pillars];
		}
		mean_i = mean(pillar_i);
		means(i % dim[0], i / dim[0]) = mean_i;
	}
	return means;
}

//' @rdname MeanPillars
// [[Rcpp::export]]
NumericMatrix VarPillars(NumericVector mat3d) {
  IntegerVector dim = mat3d.attr("dim");
	int n_pillars = dim[0] * dim[1];
	int pillar_len = dim[2];
	NumericMatrix vars(dim[0], dim[1]);
	double var_i;
	NumericVector pillar_i(pillar_len);
	for (int i = 0; i < n_pillars; i++) {
		for (int j = 0; j < pillar_len; j++) {
			pillar_i[j] = mat3d[i + j * n_pillars];
		}
		var_i = var(pillar_i);
		vars(i % dim[0], i / dim[0]) = var_i;
	}
	return vars;
}

//' @rdname MeanPillars
// [[Rcpp::export]]
NumericMatrix MedianPillars(NumericVector mat3d) {
  IntegerVector dim = mat3d.attr("dim");
	int n_pillars = dim[0] * dim[1];
	int pillar_len = dim[2];
	NumericMatrix meds(dim[0], dim[1]);
	double med_i;
	NumericVector pillar_i(pillar_len);
	for (int i = 0; i < n_pillars; i++) {
		for (int j = 0; j < pillar_len; j++) {
			pillar_i[j] = mat3d[i + j * n_pillars];
		}
		med_i = median(pillar_i);
		meds(i % dim[0], i / dim[0]) = med_i;
	}
	return meds;
}

//' Image median filter with options for dealing with NAs
//'
//' This is an alternative to \link[EBImage]{EBImage}'s
//' \code{\link[EBImage]{medianFilter}}. \code{MedianFilterB} has more options
//' for dealing with NA values. Whilst \code{\link[EBImage]{medianFilter}} can
//' either ignore NAs or set the output of any median calculation involving an
//' \code{NA} to \code{NA}, \code{MedianFilterB} can deal with \code{NA}s
//' depending on how many of them there are in a given median calculation.
//'
//' The behavior at image boundaries is such as the source image has been padded
//' with pixels whose values equal the nearest border pixel value.
//'
//' @param mat A matrix (representing an image).
//' @param size An integer; the median filter radius.
//' @param na_rm Should \code{NA}s be ignored?
//' @param na_count If this is TRUE, in each median calculation, if the majority
//' of arguments are \code{NA}s, \code{NA} is returned but if the \code{NA}s are
//' in the minority, they are ignored as in \code{median(x, na.rm = TRUE)}.
//'
//' @return A matrix (the median filtered image).
//'
//' @examples
//' m <- matrix(1:9, nrow = 3)
//' m[2:3, 2:3] <- NA
//' print(m)
//' MedianFilterB(m)
//' MedianFilterB(m, na_rm = TRUE)
//' MedianFilterB(m, na_count = TRUE)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix MedianFilterB(NumericMatrix mat, int size = 1,
                           bool na_rm = false, bool na_count = false) {
	int nr = mat.nrow();
	int nc = mat.ncol();
	NumericMatrix median_filtered(nr, nc);
	int square_side_len = 2 * size + 1;
	NumericMatrix square(square_side_len, square_side_len);
	int row;
	int col;
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < nc; j++) {
			for (int k = -size; k <= size; k++) {
				for (int l = -size; l <= size; l++) {
					row = i + k;
					col = j + l;
					/* The behavior at image boundaries is such as the source image
					has been padded with pixels whose values equal the nearest
					border pixel value. */
					while (row < 0 || row >= nr) {
						row -= boost::math::sign(row);
					}
					while (col < 0 || col >= nc) {
						col -= boost::math::sign(col);
					}
					square(k + size, l + size) = mat(row, col);
				}
			}
			double filtered_ij;
			if (na_count) {
			  if (sum(is_na(square)) > pow(square_side_len, 2) / 2.0)
			    filtered_ij = NA_REAL;
			  else
			    filtered_ij = median(square, true);
			}
			else
				filtered_ij = median(square, na_rm);
			median_filtered(i, j) = filtered_ij;
		}
	}
	return median_filtered;
}
