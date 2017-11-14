MatrixPlotB <- function(mat, orig.img = NULL, med.scale = 1,
                        scale.name = "scale", ranges = NULL, colours = NULL,
                        limits = NULL, low = "blue", high = "red",
                        clip = FALSE, clip.low = FALSE, clip.high = FALSE,
                        legend.title.size = 40, legend.label.size = 30,
                        barheight = 15) {
  plain.theme <- theme(axis.ticks = element_blank(), axis.text = element_blank(),
                       axis.title = element_blank(),
                       legend.justification = 0.5,
                       panel.background = element_rect(fill = "white"))
  if (!is.null(orig.img)) {
    orig.img.med.adj <- median(orig.img) * med.scale
    orig.img.medpillars <- MedianPillars(orig.img)
    brightness.med.adj <- median(mat) * med.scale
    to.exclude <- orig.img.medpillars <= orig.img.med.adj &
      mat <= brightness.med.adj
    mat[to.exclude] <- NA
  }
  df <- reshape2::melt(mat) %>%
    transmute(x = 1 + max(Var1) - Var1, y = Var2, value = value)
  if (!is.null(ranges)) {
    if (is.null(colours)) {
      colours <- length(ranges) - 1
    } else if (length(colours) != length(ranges) - 1) {
      stop("The number of colours must match the number of ranges")
    }
  }
  if (is.numeric(colours)) {
    if (length(colours) != 1) stop("If colours is numeric it must have length 1")
    topocols <- topo.colors(colours)
  }
  if (!is.null(colours)) {
    if (is.null(ranges)) {
      ranges <- seq(min(df$value), max(df$value),
                    length.out = length(colours) + 1)
    }
    ranges <- ranges %>% {
      cbind(.[-length(.)], .[-1])  # adjacent pairs
    }
    ranges.typeset <- apply(ranges, 1,
                            function(x) paste(round(x, 2), collapse = "-"))
    colours.ranges <- factor(df$value %>% sapply(WhichInterval, ranges)) %T>%
    {levels(.) = ranges.typeset}
    df <- df %>% mutate(colour = colours.ranges)
    ggplot(df, aes(x, y, fill = colour)) +
      scale_fill_manual(scale.name, values = topocols,
                        guide = guide_legend(
                          title.theme = element_text(size = legend.title.size, angle = 90),
                          title.position = "right",
                          title.hjust = 0.195,
                          barheight = barheight,
                          label.theme = element_text(size = legend.label.size, angle = 0))) +
      coord_fixed() +
      geom_raster() +
      plain.theme
  } else {
    if (is.null(limits)) limits <- range(df$value)
    if (clip) {
      clip.low <- TRUE
      clip.high <- TRUE
    }
    if (clip.low) df$value[df$value < limits[1]] <- limits[1]
    if (clip.high) df$value[df$value > limits[2]] <- limits[2]
    ggplot(df, aes(x, y, fill = value)) +
      scale_fill_gradient(scale.name, limits = limits,
                          low = low, high = high, na.value = "black",
                          guide = guide_colourbar(
                            title.theme = element_text(size = legend.title.size, angle = 90),
                            title.position = "right",
                            title.hjust = 0.5,
                            label.theme = element_text(size = legend.label.size, angle = 0),
                            barheight = barheight
                          )
      ) +
      coord_fixed() +
      geom_raster() +
      plain.theme
  }
}

DN <- function(img) EBImage::display(normalize(img), method = "r")

#' Subtract one side of a distribution from both sides.
#'
#' Given a (user-specified) midpoint of a distribution, subtract one side of it
#' from the other or both (the default, making one side zero) side(s), subtracting from the other side by mirroring across the midpoint. Beware that this function does not normalize the output such that it integrates to 1.
#'
#' @param samp An empirical sample of the distribution.
#' @param midpoint To say what we mean by a "side" of the distribution, we need
#'   a middle to say "the left is anything to the left of the middle and the
#'   right is anything to the right". This \code{midpoint} says where the
#'   "middle" is.
#' @param side Which side do you wish to subtract, left or right? It's good
#'   enough to use the first letter.
#' @param both Subtract the selected side from both sides of the distribution (\code{both = TRUE}), or just the other side (\code{both = FALSE})?
#' @param negs Allow negative values in the probability density of the result?
#'
#' @return A data frame with two columns x and y giving the probability density of the output. The probability density at point x is y. Beware that this is not a real probability distribution in that it doesn't integrate to 1 and indeed it allows negative probability densities.
#'
#' @examples
#' sampl <- 0:10
#' plot(density(sampl))
#' subside3 <- SubSide(sampl, 4)
#' plot(subside3, type = "l")
#' plot(SubSide(sampl, 3, both = FALSE), type = "l")
SubSide <- function(samp, midpoint, side = "left", both = TRUE, negs = FALSE) {
  side <- tolower(side)
  if (startsWith("left", side)) side <- "left"
  if (startsWith("right", side)) side <- "right"
  stopifnot(side %in% c("left", "right"))
  stopifnot(midpoint < max(samp) && midpoint > min(samp))
  dens <- stats::density(samp)
  x <- dens$x
  y <- dens$y
  midpoint <- Closest(midpoint, x)
  nleft <- sum(x < midpoint)
  nright <- sum(x > midpoint)
  difference <- nleft - nright
  zeros <- rep(0, abs(difference))
  if (difference > 0) {
    y <- c(y, zeros)
    rightmostpoint <- max(x)
    leftmostpoint <- midpoint - (rightmostpoint - midpoint)
    x <- seq(leftmostpoint, rightmostpoint, length.out = length(y))
  } else if (difference < 0) {
    y <- c(zeros, y)
    leftmostpoint <- min(x)
    rightmostpoint <- midpoint + (midpoint - leftmostpoint)
    x <- seq(leftmostpoint, rightmostpoint, length.out = length(y))
  }
  midpoint <- Closest(midpoint, x)  # avoid rounding errors below
  stopifnot(length(x) %% 2 == 1 && length(y) %% 2 == 1)
  if (side == "left") {
    lefts <- y[x < midpoint]
    y[x > midpoint] <- y[x > midpoint] - rev(lefts)
    if (both) y[x <= midpoint] <- 0
  } else {
    rights <- y[x > midpoint]
    y[x < midpoint] <- y[x < midpoint] - rev(rights)
    if (both) y[x >= midpoint] <- 0
  }
  if (!negs) y[y < 0] <- 0
  data.frame(x = x, y = y)
}

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

#' Convert a binary image to tiff.
#'
#' Read a binary image (pixels must be integers), specifying its dimension and
#' write it back to disk as a tifff file. The `Bin2TiffFolder` function does
#' this for an entire folder.
#'
#' When using `Bin2TiffFolder`, the images must all have the same dimension.
#'
#' @param bin.path The path to the binary file.
#' @param bits How many bits are there per pixel? This should be a multiple of
#'   8.
#' @param height The height of the image (in pixels).
#' @param width The width of the image (in pixels).
#' @param frames How many frames are there in the stack (default 1).
#' @param endian Eigher `"big"` or `"little"` (see [readBin]).
#' @param signed Logical. Only used for integers of sizes 1 and 2, when it
#'   determines if the quantity on file should be regarded as a signed or
#'   unsigned integer.
#'
#' @export
Bin2Tiff <- function(bin.path, bits, height, width, frames = 1,
                     endian = .Platform$endian, signed = TRUE) {
  path.short <- bin.path
  if (stringr::str_detect(path.short, "/")) {
    file.dir <- filesstrings::str_before_nth(path.short, "/", -1)
    current.wd <- getwd()
    on.exit(setwd(current.wd))
    setwd(file.dir)
    path.short <- stringr::str_split(path.short, "/") %>%
      unlist %>%
      dplyr::last()
  }
  bin.vector <- readBin(path.short, "int", n = height * width * frames,
                        signed = signed, size = bits %/% 8)
  arr <- array(bin.vector, dim = c(height, width, frames))
  if (dim(arr)[3] == 1) arr <- arr[, , 1]
  WriteIntImage(arr, filesstrings::before_last_dot(path.short))
}

#' @rdname Bin2Tiff
#'
#' @param folder.path The path to the folder (defaults to current location).
#'
#' @examples
#' setwd(tempdir())
#' dir.create("temp_dir")
#' file.copy(system.file("extdata", "b.bin", package = "nandb"), "temp_dir")
#' Bin2Tiff("temp_dir/b.bin", height = 2, width = 2, bits = 8)
#' Bin2TiffFolder("temp_dir", height = 2, width = 2, bits = 8)
#' list.files("temp_dir")
#' @export
Bin2TiffFolder <- function(folder.path = ".", bits, height, width, frames = 1,
                           endian = .Platform$endian, signed = TRUE) {
  current.wd <- getwd()
  on.exit(setwd(current.wd))
  setwd(folder.path)
  bin.names <- list.files(pattern = "\\.bin$")
  lapply(bin.names, Bin2Tiff, bits, height, width, frames, endian, signed)
}

CanBeInteger <- function(x) floor(x) == x

WhichInterval <- function(numbers, interval.mat) {
  # interval.mat is a two-column matrix where the rows are
  # increasing, non-intersecting, half-open (open at the right)
  # intervals on the real line
  sq <- as.vector(t(interval.mat))  # transposed interval matrix
  if (!(all(diff(sq) >= 0) && all(interval.mat[, 1] < interval.mat[, 2]))) {
    stop("interval.mat must be a two-column matrix where the rows are ",
         "increasing, non-intersecting, half-open intervals on the real line.")
  }
  WhichIntervalC(numbers, interval.mat)
}

#' Turn a 3d array into a list of pillars
#'
#' Suppose the array is of dimension `n1 * n2 * n3`, then pillar \code{[i,
#' j, ]} is list element `i + n1 * (j - 1)`, so the first element is pillar
#' `[1, 1, ]`, the second is pillar `[1, 2, ]` and so on.
#'
#' @param arr3d A 3-dimensional array.
#'
#' @return A list.
#'
#' @seealso [PillarsListToArr]
#'
#' @examples
#' arr <- array(1:27, dim = rep(3, 3))
#' print(arr)
#' ListPillars(arr)
#'
#' @export
ListPillars <- function(arr3d) {
  d <- dim(arr3d)
  lapply(seq_len(prod(d[1:2])), function(x) {
    arr3d[((x - 1) %% d[1]) + 1, ((x - 1) %/% d[1]) + 1, ]
  })
}

#' Make a list of pillars back into a 3D array.
#'
#' #' Suppose the array is of dimension `n1 * n2 * n3`, then pillar
#' `[i, j, ]` is list element `i + n1 * (j - 1)`, so the first element
#' is pillar `[1, 1, ]`, the second is pillar `[1, 2, ]` and so on.
#'
#' @param pillars.list The list of pillars
#' @param dim The desired dimension of the output array.
#'
#' @return A 3-dimensional array.
#'
#' @seealso [ListPillars]
#'
#' @examples
#' arr <- array(1:27, dim = rep(3, 3))
#' print(arr)
#' ListPillars(arr)
#' PillarsListToArr(ListPillars(arr), dim(arr))
#'
#' @export
PillarsListToArr <- function(pillars.list, dim) {
  arr3d <- array(0, dim = dim)
  for (i in seq_along(pillars.list)) {
    arr3d[((i - 1) %% dim[1]) + 1,
          ((i - 1) %/% dim[1]) + 1,
          ] <- pillars.list[[i]]
  }
  arr3d
}

#' Make each pillar of a 3D array into a column of a tibble
#'
#' Create a data frame (tibble) where pillar `[i, j, ]` is column \code{i + n1 *
#' (j - 1)} with column name "i_j" (we use an underscore here rather than a
#' comma because a comma is the delimiter in csv files so writing this data
#' frame to a csv could cause confusion).
#'
#' @param arr3d A 3-dimensional array.
#'
#' @return A [tibble][tibble::tibble].
#'
#' @examples
#' arr <- array(1:27, dim = rep(3, 3))
#' print(arr)
#' ListPillars(arr)
#' PillarsDF(arr)
#'
#' @export
PillarsDF <- function(arr3d) {
  pillars.list <- ListPillars(arr3d)
  d <- dim(arr3d)
  namez <- expand.grid(seq_len(d[1]), seq_len(d[2])) %>% {
    stringr::str_c(.$Var1, "_", .$Var2)
  }
  names(pillars.list) <- namez
  tibble::as_tibble(pillars.list)
}

BrightnessVec <- function(vec) {
  if (filesstrings::all_equal(vec, 0)) return(0)
  stats::var(vec)/mean(vec)
}

#' Fix lookup table error when reading images.
#'
#' Even if an image stored on disk has 1 channel only, if it has an associated
#' LUT (lookup table, to tell programs like ImageJ to display them in (for
#' example) green rather than grey, then [EBImage::readImage()] and
#' hence [ReadImageData()] will read in the image as a 3-channel rgb
#' colour image, where two of the channels have all-zero intensity values and
#' one of them has the values you wanted. This function deletes the zero
#' channels from the array. If there were no such LUT errors and the image read
#' in as you desired, then this function does nothing.
#'
#' @param arr An array, representing the read image.
#' @param ndim.out How many dimensions do you want the output array to have?
#'
#' @return An array. If the function implemented a 'fix', then the output array
#'   will have one less dimension than the input array.
#'
#' @examples
#' has.lut.error <- abind::abind(matrix(0, nrow = 2, ncol = 2),
#'                               matrix(1:4, nrow = 2),
#'                               matrix(0, nrow = 2, ncol = 2),
#'                               along = 3)
#' has.lut.error
#' FixLUTError(has.lut.error, 2)
#' has.lut.error3d <- abind::abind(has.lut.error, has.lut.error, along = 4)
#' FixLUTError(has.lut.error3d, 3)
#' @export
FixLUTError <- function(arr, ndim.out) {
  d <- dim(arr)
  nd <- length(d)
  if (nd == ndim.out)
    return(arr)
  error.msg <- paste("Modification is needed to get the output array to have",
                     "the output array have the number of dimensions that you require,",
                     "however this function does not know how to make those modifications.")
  if (nd - 1 != ndim.out) {
    error.msg <- paste0("This function can only modify your image to reduce ",
                        "its dimensionality by 1, however R is reading in your image as having ",
                        nd, " dimensions, whereas you wish your image to be in a form whereby it",
                        " has ", ndim.out, " dimensions, this is a dimension reduction of ",
                        nd - ndim.out, ".")
    stop(error.msg)
  }
  if (d[3] != 3) {
    err.msg <- paste0("This function expects that in the read image, ",
                      "the third array slot is the rgb colour slot and hence ",
                      "that the third dimension has value 3, however in the ",
                      "read image, this slot had dimension value ",
                      d[3], ".")
    stop(err.msg)
  }
  third.solos <- lapply(seq_len(d[3]), function(x) {
    R.utils::extract(arr, seq_len(d[1]), seq_len(d[2]), x,
                     drop = TRUE)
  })
  nonzero <- !vapply(lapply(third.solos, as.vector), filesstrings::all_equal,
                     logical(1), 0)
  if (sum(nonzero) != 1) {
    err.msg <- paste0("This function expects that in the read image, ",
                      "the third array slot is the rgb colour slot and that ",
                      "all but one of these color channels is zero, ",
                      "however in the read image, ",
                      sum(nonzero), " of these were nonzero.")
    stop(err.msg)
  }
  R.utils::extract(arr, seq_len(d[1]), seq_len(d[2]), nonzero,
                   drop = TRUE)
}


#' What's the closest value in a vector?
#'
#' Given a number `x` and a numeric vector `vec`, what's the closest value to
#' `x` in `vec`?
#'
#' @param x A number.
#' @param vec A numeric vector.
#' @param index If set to `TRUE`, return the index (rather than the value) of
#'   the closest element.
#'
#' @return A number.
#'
#' @examples
#' Closest(pi, 0:10)
#' Closest(pi, 0:10, index = TRUE)
#'
#' @export
Closest <- function(x, vec, index = FALSE) {
  ind <- which.min(abs(x - vec))
  ifelse(index, ind, vec[ind])
}

#' Collapse a big set of ranges into a smaller set.
#'
#' Say you have many ranges (or bins) in which you're assigning continuous
#' values and you'd like to collapse these ranges such that there are fewer of
#' them (but they still cover the same part of that continuous scale). This
#' function is here to help.
#'
#' One property of this procedure is that each new range is the union of old
#' ranges.
#'
#' @param ranges A strictly increasing (numeric) vector. Each set of adjacent
#'   elements are interpreted as the bounds of a range (or bin).
#' @param n.out A natural number. How many ranges should the output have?
#' @param preserve A vector. Are there any original ranges that you'd like to
#'   preserve? If so set them here. The first range is the interval
#'   \code{[ranges[1], ranges[2])} and so on.
#' @param prefer.low Are you more interested in the lower ranges? If so, set
#'   this to true and all the high ranges will be collapsed into one.
#' @param prefer.high Are you more interested in the higher ranges? If so, set
#'   this to true and all the high ranges will be collapsed into one.
#'
#' @return A vector of the new ranges.
#'
#' @examples
#' set.seed(0)
#' ranges <- sort(sample(1:100, 11))
#' print(ranges)
#' CollapseRanges(ranges, 6, c(3, 7))
#' CollapseRanges(ranges, 6, c(3, 7), prefer.low = TRUE)
#' CollapseRanges(ranges, 6, c(3, 7), prefer.high = TRUE)
#' @export
CollapseRanges <- function(ranges, n.out, preserve = NULL, prefer.low = FALSE,
                           prefer.high = FALSE) {
  stopifnot(is.vector(ranges), is.numeric(ranges))
  if (prefer.high && prefer.low) {
    stop("One cannot select both prefer.high and prefer.low.")
  }
  if (!all(diff(ranges) > 0)) stop("ranges must be strictly increasing.")
  if (n.out >= length(ranges)) {
    stop("Cannot collapse ", length(ranges) - 1, " ranges into ", n.out,
         " ranges. Collapsing must reduce the number of ranges.")
  }
  lp <- length(preserve)
  if (lp > n.out) {
    stop("One cannot try to preserve more ranges than one wants overall. ",
         "i.e. one cannot have length(preserve) > n.out.")
  }
  between.preserve <- GroupClose(setdiff(seq_along(ranges),
                                         preserve))
  lbp <- length(between.preserve)
  nrbp <- length(unlist(between.preserve)) - 1
  ranges.lower.bound <- lp + lbp
  if (ranges.lower.bound > n.out) {
    message <- paste0("The way in which you've chosen n.out and preserve makes",
                      " it impossible to collapse the entirety of ranges into ",
                      "n.out = ", n.out, " ranges. your selection for preserve",
                      " leaves ", lp, " ranges to preserve and ", lbp,
                      " sets of ranges in between, so the result must have at ",
                      "least ", ranges.lower.bound, " ranges, which is more ",
                      "than you want.")
    stop(message)
  }
  required.nrbp <- n.out - lp
  req.nrbp.reduction <- nrbp - required.nrbp
  lsbp <- lengths(between.preserve)
  to.merge <- between.preserve %T>% {
    .[[length(.)]] <- .[[length(.)]][-length(.[[length(.)]])]
  } %>% lapply(function(x) as.list(seq_along(x)))
  lstm <- lengths(to.merge)
  if (prefer.high) {
    max.reductions <- lstm - 1
    reductions.to.make <- 0
    i <- 1
    while (req.nrbp.reduction > 0) {
      amount.i <- min(req.nrbp.reduction, max.reductions[i])
      reductions.to.make[i] <- amount.i
      req.nrbp.reduction <- req.nrbp.reduction - amount.i
      i <- i + 1
    }
    for (i in seq_along(reductions.to.make)) {
      rtm.i <- reductions.to.make[i]
      if (rtm.i > 0) {
        slrtm.ip1 <- seq_len(rtm.i + 1)
        to.merge[[i]][[1]] <- unlist(to.merge[[i]][slrtm.ip1])
        to.merge[[i]][slrtm.ip1[-1]] <- NULL
      }
    }
  } else if (prefer.low) {
    max.reductions <- lstm - 1
    reductions.to.make <- rep(0, length(max.reductions))
    i <- length(max.reductions)
    while (req.nrbp.reduction > 0) {
      amount.i <- min(req.nrbp.reduction, max.reductions[i])
      reductions.to.make[i] <- amount.i
      req.nrbp.reduction <- req.nrbp.reduction - amount.i
      i <- i - 1
    }
    for (j in seq(length(reductions.to.make), i + 1)) {
      rtm.j <- reductions.to.make[j]
      if (rtm.j > 0) {
        ltmj <- length(to.merge[[j]])
        start.merge.j <- ltmj - rtm.j
        sq <- seq(start.merge.j, ltmj)
        to.merge[[j]][[start.merge.j]] <- unlist(to.merge[[j]][sq])
        to.merge[[j]][sq[-1]] <- NULL
      }
    }
  } else {
    between.preserve.range.lengths <- lapply(between.preserve,
                                             function(x) {
                                               ranges[x + 1] - ranges[x]
                                             }) %>% lapply(stats::na.omit)
    bpls.new <- between.preserve.range.lengths
    while (sum(lengths(bpls.new)) > required.nrbp) {
      combined.adjacent.lengths <- lapply(bpls.new, RcppRoll::roll_sum,
                                          2) %T>% {
                                            .[lengths(.) < 1] <- NA
                                          }
      min.group <- which.min(vapply(combined.adjacent.lengths, min,
                                    numeric(1)))
      min.ingroup.index <- which.min(combined.adjacent.lengths[[min.group]])
      to.merge[[min.group]][[min.ingroup.index]] <-
        c(to.merge[[min.group]][[min.ingroup.index]],
          to.merge[[min.group]][[min.ingroup.index + 1]])
      to.merge[[min.group]][[min.ingroup.index + 1]] <- NULL
      bpls.new[[min.group]][min.ingroup.index] <-
        min(unlist(combined.adjacent.lengths), na.rm = TRUE)
      bpls.new[[min.group]] <- bpls.new[[min.group]][-(min.ingroup.index + 1)]
    }
  }
  merging.range.nums <- mapply(function(x, y) {
    lapply(y, function(a) x[a])
  }, between.preserve, to.merge, SIMPLIFY = FALSE) %>% c(as.list(preserve)) %>%
    Reduce(c, .) %>% {
      .[order(vapply(., min, numeric(1)))]
    }
  merging.range.nums %>% lapply(function(x) ranges[seq(x[1],
                                                       dplyr::last(x) + 1)]) %>% lapply(range) %>% unlist %>%
    unique
}

#' Group together close adjacent elements of a vector.
#'
#' Given a strictly increasing vector (each element is bigger than the last),
#' group together stretches of the vector where \emph{adjacent} elements are
#' separeted by at most some specified distance. Hence, each element in each
#' group has at least one other element in that group that is \emph{close} to
#' it. See the examples.
#' @param vec.ascending A strictly increasing numeric vector.
#' @param max.gap The biggest allowable gap between adjacent elements for them
#'   to be considered part of the same \emph{group}.
#' @return A where each element is one group, as a numeric vector.
#' @examples
#' GroupClose(1:10, 1)
#' GroupClose(1:10, 0.5)
#' GroupClose(c(1, 2, 4, 10, 11, 14, 20, 25, 27), 3)
#' @export
GroupClose <- function(vec.ascending, max.gap = 1) {
  lv <- length(vec.ascending)
  if (lv == 0)
    stop("vec.ascending must have length greater than zero.")
  test <- all(vec.ascending > dplyr::lag(vec.ascending), na.rm = TRUE)
  if (!test)
    stop("vec.ascending must be strictly increasing.")
  if (lv == 1) {
    return(list(vec.ascending))
  } else {
    gaps <- vec.ascending[2:lv] - vec.ascending[1:(lv - 1)]
    big.gaps <- gaps > max.gap
    nbgaps <- sum(big.gaps)  # number of big (>10) gaps
    if (!nbgaps) {
      return(list(vec.ascending))
    } else {
      ends <- which(big.gaps)  # vertical end of lines
      group1 <- vec.ascending[1:ends[1]]
      lg <- list(group1)
      if (nbgaps == 1) {
        lg[[2]] <- vec.ascending[(ends[1] + 1):lv]
      } else {
        for (i in 2:nbgaps) {
          lg[[i]] <- vec.ascending[(ends[i - 1] + 1):ends[i]]
          ikeep <- i
        }
        lg[[ikeep + 1]] <- vec.ascending[(ends[nbgaps] +
                                            1):lv]
      }
      return(lg)
    }
  }
}

#' Get the best \emph{spread} of numbers in an interval.
#'
#' Say we have an interval \eqn{[a, b]}  and we want to evenly spread \eqn{n}
#' numbers in this interval, however \eqn{m} of these \eqn{n} numbers have been
#' chosen beforehand, so we have to fithe the other \eqn{n - m} in around them
#' the best we can. Our goal is to maximize the minimum distance between
#' adjacent elements.
#'
#' The idea is to, for an interval of size \eqn{s}, put in \eqn{floor(s / (n -
#' m))} numbers into each of these intervals. Then, if there any numbers left
#' over, put exactly one into each of the intervals with the biggest \eqn{(s /
#' (n - m)) - `floor`(s / (n - m))} (starting with the biggest and working
#' one's way down) until there are no numbers left to assign. The end intervals
#' need special treatment since (assuming the ends are not in the specific
#' numbers), the end intervals are open on one side (one boundary has no point
#' on it), whereas all the middle intervals are not.
#'
#' @param interval A length 2 numeric vector. The real interval over which one
#'   wants to spread the numbers.
#' @param specific A numeric vector. The specific numbers that one wants to be
#'   part of our spread.
#' @param n A number. The total number of numbers that you want to be in the
#'   spread (including the number(s) in `specific`).
#' @param log Measure the difference between adjacent elements as the difference
#'   in their log instead?
#'
#' @return A numeric vector. The chosen \eqn{n} numbers.
#'
#' @examples
#' SpreadSpecific(c(0, 10), 1, 3)
#'
#' @export
SpreadSpecific <- function(interval, specific, n, log = FALSE) {
  stopifnot(length(interval) == 2, length(specific) > 0, length(n) == 1)
  interval <- sort(interval)
  if (log) {
    if (interval[1] <= 0) {
      stop("If log is selected, ",
           "the interval must be on the positive real line")
    }
  }
  if (any(specific <= interval[1]) | any(specific >= interval[2])) {
    stop("All members of specific must fall in interval.")
  }
  specific <- sort(unique(specific))
  lspec <- length(specific)
  interval.pops.init <- c(1, rep(2, lspec - 1), 1)
  intervals <- unique(c(interval[1], specific, interval[2]))
  if (log) {
    interval.lengths <- diff(log(intervals))
  } else {
    interval.lengths <- diff(intervals)
  }
  stopifnot(length(interval.lengths) == length(interval.pops.init))
  interval.pops.final <- SpreadSpecificHelper(interval.lengths,
                                              interval.pops.init, n - lspec)
  interval.pops.add <- interval.pops.final - interval.pops.init
  intervals.list <- cbind(dplyr::lag(intervals), intervals)[-1, ] %>%
    Mat2RowList
  to.be.added.to <- interval.pops.add > 0
  intervals.list.to.add <- intervals.list[to.be.added.to]
  interval.pops.final.to.add <- interval.pops.final[to.be.added.to]
  if (log) {
    for (i in seq_along(intervals.list.to.add)) {
      intervals.list.to.add[[i]] <- exp(seq(log(intervals.list.to.add[[i]][1]),
                                            log(intervals.list.to.add[[i]][2]),
                                            length.out = interval.pops.final.to.add[i]))
    }
  } else {
    for (i in seq_along(intervals.list.to.add)) {
      intervals.list.to.add[[i]] <- seq(intervals.list.to.add[[i]][1],
                                        intervals.list.to.add[[i]][2],
                                        length.out = interval.pops.final.to.add[i])
    }
  }
  out <- sort(unique(c(specific, unlist(intervals.list.to.add))))
  diff.out <- diff(out)
  below.tol <- diff.out < 1.5e-8
  out[!below.tol]
}

ScaleBreaksFun <- function(include.breaks = NULL, log = FALSE) {
  if (is.null(include.breaks) && (!log)) return(waiver())
  function(x) {
    if (is.null(include.breaks)) {
      if (log) {
        untruncated <- exp(seq(log(min(x)), log(max(x)), length.out = 5))
      }
    } else {
      untruncated <- SpreadSpecific(x, include.breaks, 5, log = log)
    }
    min.sep <- min(diff(untruncated))
    roundn <- log10(min.sep) %>% {
      # get a sensible amount of digits for breaks
      ifelse(. >= 1, 0, -floor(.))
    }
    brks <- ifelse(untruncated %in% include.breaks,
                   untruncated, round(untruncated, roundn))
    # if (brks[1] < min(x)) brks[1] <- brks[1] + 10 ^ -roundn
    lbrks <- length(brks)
    # if (brks[lbrks] > max(x)) brks[lbrks] <- brks[lbrks] - 10 ^ -roundn
    brks
  }
}

bpp <- function(mcc, seed = NULL) {
  suppressWarnings(BiocParallel::MulticoreParam(workers = mcc, RNGseed = seed))
}

ListChannels <- function(arr4d, n.ch) {
  d <- dim(arr4d)
  stopifnot(length(d) == 4)
  if (d[3] != n.ch) {
    stop("You have input an image with ", d[3], " channels, ",
         "but you have indicated that it has ", n.ch, " channels.")
  }
  purrr::map(seq_len(d[3]), ~ arr4d[, , ., ])
}

ChannelList2Arr <- function(chnl.lst) {
  stopifnot(is.list(chnl.lst))
  dims <- purrr::map(chnl.lst, dim)
  if (!filesstrings::all_equal(dims)) {
    stop("The dimensions of the elements of chnl.lst must all be the same.")
  }
  if (length(dims[[1]]) == 2) {
    out <- purrr::reduce(chnl.lst, ~ abind::abind(.x, .y, along = 3))
  } else if (length(dims[[1]]) == 3) {
    out <- purrr::reduce(chnl.lst, ~ abind::abind(.x, .y, along = 4)) %>%
      aperm(c(1, 2, 4, 3))
  } else {
    stop("The elements of chnl.lst must be 2- or 3-dimensional.")
  }
  if (all(purrr::map_lgl(chnl.lst, ~ "frames" %in% names(attributes(.))))) {
    frames <- purrr::map_int(chnl.lst, ~ attr(., "frames")) %>% unique()
    if (length(frames) == 1) attr(out, "frames") <- frames
  }
  att_names <- purrr::map(chnl.lst, ~ names(attributes(.)))
  for (att_name in c("tau", "mst", "filt")) {
    if (all(purrr::map_lgl(att_names, ~ att_name %in% .))) {
      atts <- purrr::map(chnl.lst, ~ attr(., att_name))
      if (filesstrings::all_equal(atts)) {
        attr(out, att_name) <- paste(unlist(atts), collapse = ",")
      }
    }
  }
  out
}

## Translate the fail argument from what the user selects to what failed pixels
## should be set to.
TranslateFail <- function(arr, fail) {
  stopifnot(length(fail) == 1)
  if (is.na(fail)) return(NA)
  if (is.numeric(fail)) {
    stopifnot(fail >= 0)
  } else if (is.character(fail)) {
    fail <- RSAGA::match.arg.ext(fail, c("saturate", "zero"),
                                 ignore.case = TRUE)
  }
  if (fail == "zero") {
    fail = 0
  } else if (fail == "saturate") {
    mx <- max(arr, na.rm = TRUE)
    if (mx >= 2 ^ 8) {
      bits.per.sample <- 16
    } else {
      bits.per.sample <- 8
    }
    fail <- 2 ^ bits.per.sample - 1
  }
  if (is.numeric(fail) && !is.na(fail)) fail <- as.integer(round(fail))
  fail
}
