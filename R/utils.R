CanBeInteger <- function(x) floor(x) == x

Slices <- function(slices, mat3d) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (!all(slices %in% seq_len(d[3]))) stop("slices out of range of slices of mat3d")
  mat3d[, , slices]
}
WhichInterval <- function(numbers, interval.mat) {
  # interval.mat is a two-column matrix where the rows are increasing, non-intersecting, half-open (open at the right) intervals on the real line
  tim <- t(interval.mat)  # transposed interval matrix
  if (!all(tim[-length(tim)] <= tim[-1]) && all(interval.mat[, 1] < interval.mat[, 2])) {
    stop("interval.mat must be a two-column matrix where the rows are increasing, non-intersecting, half-open intervals on the real line")
  }
  WhichIntervalC(numbers, interval.mat)
}
AllEqual <- function(a, b = NA, allow = T, cn = F) {
  if (is.na(b[1])) {
    return(length(unique(a)) == 1)
  } else {
    if (allow) {
      if (length(a) == 1) {
        a <- rep(a, length(b))
        if (is.array(b)) b <- as.vector(b)
      }
      if (length(b) == 1) {
        b <- rep(b, length(a))
        if (is.array(a)) a <- as.vector(a)
      }
    }
    return(isTRUE(all.equal(a, b, check.names = cn)))
  }
}

#' Apply a function to each pillar of a 3-dimensional array.
#'
#' Define a "pillar" of a 3-dimensional array as pillar \code{i,j} off array \code{arr} being \code{arr[i, j, ]}. This function applies a specified function to each pillar.
#'
#' @param mat3d A 3-dimensional array.
#' @param FUN A function which takes a vector as imput and, for a given input length, outputs a vector of constant length (can be 1).
#'
#' @return If \code{FUN} is returning length 1 vectors, a matrix whereby \code{mat[i, j] = FUN(mat3d[i, j, ])}. If FUN is returning vectors of length \code{l > 1}, a 3-dimensional array whereby \code{arr[i, j, ] = FUN(mat3d[i, j, ])}.
#' @export
ApplyOnPillars <- function(mat3d, FUN) {
  apply(mat3d, c(1, 2), FUN) %>%  # This applies FUN over pillars of mat3d (c(1, 2) means iterate over the first two slots of mat3d so FUN is recursively applied to mat3d[i, j, ] (notice i and j occupying the first two slots there) for all possible i, j) where the result from pillar [i, j, ] goes into column [, i, j] of the result.
    (function(x) {
      if (length(dim(x)) == 3) {
        return(aperm(x, c(2, 3, 1)))  # After the apply step, we have the column [, i, j] containing what we would like to be pillar [i, j, ]. Hence we need aperm to move the second index to the first, the third to the second and the first to the third, hence the c(2, 3, 1)
      } else {
        return(x)  # If fun returns vectors of length 1, apply will return a matrix (2D) so we don't have to work any magic (it's already as we would like it)
      }
    })
}
BrightnessVec <- function(vec) {
  if (AllEqual(vec, 0)) return(0)
  stats::var(vec) / mean(vec)
}

#' Fix lookup table error when reading images.
#'
#' Even if an image stored on disk has 1 channel only, if it has an associated
#' LUT (lookup table, to tell programs like ImageJ to display them in (for
#' example) green rather than grey, then \code{\link[EBImage]{readImage}} and
#' hence \code{\link{ReadImageData}} will read in the image as a 3-channel rgb
#' colour image, where two of the channels have all-zero intensity values and
#' one of them has the values you wanted. This function deletes the zero
#' channels from the array. If there were no such LUT errors and the image read
#' in as you desired, then this function does nothing.
#'
#' @param arr An array, representing the read image.
#' @param ndim.out How many dimensions do you want the output array to have?
#'
#' @return An array. If the function implemented a "fix", then the output array
#'   will have one less dimension than the input array.
#' @export
FixLUTError <- function(arr, ndim.out) {
  d <- dim(arr)
  nd <- length(d)
  if (nd == ndim.out) return(arr)
  error.msg <- "Modification is needed to get the output array to have the output array have the number of dimensions that you require, however this function does not know how to make those modifications."
  if (nd - 1 != ndim.out) {
    error.msg <- paste0("This function can only modify your image to reduce its dimensionality by 1, however R is reading in your image as having ",
                        nd,
                        " dimensions, whereas you wish your image to be in a form whereby it has ",
                        ndim.out,
                        " dimensions, this is a dimension reduction of ",
                        nd - ndim.out,
                        ".")
    stop(error.msg)
  }
  if (d[3] != 3) {
    err.msg <- paste0("This function expects that in the read image, the third array slot is the rgb colour slot and hence that the third dimension has value 3, however in the read image, this slot had dimension value ",
                      d[3], ".")
    stop(err.msg)
  }
  third.solos <- lapply(seq_len(d[3]), function(x) {
    R.utils::extract(arr, seq_len(d[1]), seq_len(d[2]), x, drop = TRUE)
  })
  nonzero <- !sapply(third.solos, AllEqual, 0)
  if (sum(nonzero) != 1) {
    err.msg <- paste0("This function expects that in the read image, the third array slot is the rgb colour slot and that all but one of these color channels is zero, however in the read image, ", sum(nonzero), " of these were nonzero.")
    stop(err.msg)
  }
  R.utils::extract(arr, seq_len(d[1]), seq_len(d[2]), nonzero, drop = TRUE)
}

Closest <- function(x, vec) {
  which.min(abs(x - vec)) %>% {vec[.]}
}

MCLapply <- function(x, fun, ..., mcc = parallel::detectCores()) {
  dots <- list(...)
  if (tolower(.Platform$OS.type) == "windows") mcc <- 1
  if (length(dots)) {
    args <- list(X = x, FUN = fun, ... = ..., mc.cores = mcc)
  } else {
    args <- list(X = x, FUN = fun, mc.cores = mcc)
  }
  do.call(parallel::mclapply, args)
}
MCMapply <- function(fun, ..., mcc = parallel::detectCores()) {
  dots <- list(...)
  x <- dots[[1]]
  if (tolower(.Platform$OS.type) == "windows") mcc <- 1
  args <- list(FUN = fun, ... = ..., mc.cores = mcc)
  do.call(parallel::mcmapply, args)
}
