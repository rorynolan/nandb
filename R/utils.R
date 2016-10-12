CanBeInteger <- function(x) floor(x) == x

FixImage4to3 <- function(mat4d, must.start.4d = FALSE) {
  d <- dim(mat4d)
  if (must.start.4d) if (length(d) == 3) stop("mat4d is not 4-dimensional")
  if (length(d) < 4) return(mat4d)
  if (length(d) != 4) stop("mat4d must be a 4-dimensional array")
  color.sums <- rep(FALSE, d[3])
  for (i in seq_len(d[3])) color.sums[i] <- sum(mat4d[, , i, ])
  g0 <- color.sums > 0
  if (sum(g0) != 1) stop("Fixing has not worked since the input was not of the type that could be fixed")
  third.index.to.keep <- match(TRUE, g0)
  mat4d[, , third.index.to.keep, ]
}

Slices <- function(slices, mat3d) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  if (!all(slices %in% seq_len(d[3]))) stop("slices out of range of slices of mat3d")
  mat3d[, , slices]
}
WhichInterval <- function(number, interval.mat) {
  # interval.mat is a two-column matrix where the rows are increasing, non-intersecting, half-open (open at the right) intervals on the real line
  tim <- t(interval.mat)  # transposed interval matrix
  if (!all(tim[-length(tim)] <= tim[-1]) && all(interval.mat[, 1] < interval.mat[, 2])) {
    stop("interval.mat must be a two-column matrix where the rows are increasing, non-intersecting, half-open intervals on the real line")
  }
  interval <- (number >= interval.mat[, 1] & number <= interval.mat[, 2]) %>% match(T, .)
}
AllEqual <- function(a, b = NA, allow = T, cn = F) {
  if (is.na(b[1])) {
    return(length(unique(a)) == 1)
  } else {
    if (allow) {
      if (length(a) == 1) a <- rep(a, length(b))
      if (length(b) == 1) b <- rep(b, length(a))
    }
    return(isTRUE(all.equal(a, b, check.names = cn)))
  }
}
ApplyOnPillars <- function(mat3d, FUN) {
  # FUN is a function which must take a vector as an argument and for a given input length, it must have a constant output length (i.e. length(x) = length(y) => length(FUN(x)) = length(FUN(y)))
  # Say d = dim(mat3d). We wish to act on every pillar (where by pillar ij we mean mat3d[i, j, ]) and return a 3d array out.arr where all.equal(dim(out.arr)[1:2], d[1:2]) where out.arr[i, j, ] = FUN(mat3d[i, j, ])
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
  var(vec) / mean(vec)
}
MCLapply <- function(x, fun, ..., mcc = parallel::detectCores()) {
  dots <- list(...)
  if (.Platform$OS.type == "windows") mcc <- 1
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
  if (.Platform$OS.type == "windows") mcc <- 1
  args <- list(FUN = fun, ... = ..., mc.cores = mcc)
  do.call(parallel::mcmapply, args)
}
MCSapply <- function(x, fun, ..., mcc = parallel::detectCores()) {
  simplify2array(MCLapply(x, fun, ..., mcc = mcc))
}
PillarIJ <- function(mat3d, i.or.ij, j = NA) {
  if (!length(i.or.ij) %in% 1:2) stop("i.or.ij must be of length 1 or 2")
  if (is.na(j)) {
    if (length(i.or.ij) != 2) stop("if j is not specified, i.or.ij must be of length 2")
    j <- i.or.ij[2]
  } else {
    if (length(i.or.ij) != 1) stop("if j is specified, i.or.ij must be of length 1")
  }
  mat3d[i.or.ij[1], j, ]
}
PillarsList <- function(mat3d) {
  d <- dim(mat3d)
  if (length(d) != 3) stop("mat3d must be a three-dimensional array")
  ijs <- expand.grid(1:d[1], 1:d[2])
  pillars <- t(apply(ijs, 1, PillarIJ, mat3d = mat3d))
  list(pillars = pillars, coords = ijs)
}
PillarsListToMat3d <- function(pillars) {
  if (!all(names(pillars) == c("pillars", "coords"))) stop("pillars must be a list like that returned by Pillars")
  nr <- nrow(pillars$pillars)
  if (nr != nrow(pillars$coords)) stop("The elements of the pillars list must have equal numbers of rows")
  mat3d <- array(0, dim = c(apply(pillars$coords, 2, max), ncol(pillars$pillars)))  # initialise the matrix which will eventually contain the correct values
  for (i in seq_len(nr)) mat3d[pillars$coords[i, 1], pillars$coords[i, 2], ] <- pillars$pillars[i, ]
  mat3d
}
RowAdd <- function(mat, vec) {
  # add the number vec[i] to the ith row of mat for all row numbers i
  nr <- nrow(mat)
  if (length(vec) != nr) stop("The length of vec must be equal to the number of rows in mat.")
  to.add <- matrix(rep(vec, ncol(mat)), nrow = nr)
  mat + to.add
}

MyMedianFilter <- function(img.mat, size) {
  m <- max(img.mat)
  img.mat <- img.mat / m
  EBImage::medianFilter(img.mat, size) * m
}


