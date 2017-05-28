CanBeInteger <- function(x) floor(x) == x

Slices <- function(slices, mat3d) {
  d <- dim(mat3d)
  if (length(d) != 3)
    stop("mat3d must be a three-dimensional array")
  if (!all(slices %in% seq_len(d[3])))
    stop("slices out of range of slices of mat3d")
  mat3d[, , slices]
}
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

#' Apply a function to each pillar of a 3-dimensional array.
#'
#' Define a 'pillar' of a 3-dimensional array as pillar `i,j` off array
#' `arr` being `arr[i, j, ]`. This function applies a specified
#' function to each pillar.
#'
#' @param mat3d A 3-dimensional array.
#' @param FUN A function which takes a vector as imput and, for a given input
#'   length, outputs a vector of constant length (can be 1).
#'
#' @return If `FUN` is returning length 1 vectors, a matrix whereby
#'   `mat[i, j] = FUN(mat3d[i, j, ])`. If FUN is returning vectors of
#'   length `l > 1`, a 3-dimensional array whereby \code{arr[i, j, ] =
#'   FUN(mat3d[i, j, ])}.
#' @export
ApplyOnPillars <- function(mat3d, FUN) {
  apply(mat3d, c(1, 2), FUN) %>% {
    if (length(dim(.)) == 3) {
      return(aperm(., c(2, 3, 1)))
      # After the apply step, we have the column [, i, j] containing what we
      # would like to be pillar [i, j, ]. Hence we need aperm to move the second
      # index to the first, the third to the second and the first to the third,
      # hence the c(2, 3, 1)
    } else {
      return(.)  # If fun returns vectors of length 1, apply will return a
      # matrix (2D) so we don't have to work any magic
      # (it's already as we would like it)
    }
  }
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
  if (filesstrings::AllEqual(vec, 0)) return(0)
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
  nonzero <- !vapply(lapply(third.solos, as.vector), filesstrings::AllEqual,
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

Mat2ColList <- function(mat) {
  lapply(seq_len(ncol(mat)), function(i) mat[, i])
}
Mat2RowList <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) mat[i, ])
}

bpp <- function(mcc, seed = NULL) {
  suppressWarnings(BiocParallel::MulticoreParam(workers = mcc, RNGseed = seed))
}

ListChannels <- function(arr4d) {
  d <- dim(arr4d)
  stopifnot(length(d) == 4)
  if (d[3] > 4) {
    stop("You have input an image with ", d[3], "channels. ",
         "This function can deal with images of at most 4 channels.")
  }
  purrr::map(seq_len(d[3]), ~ arr4d[, , ., ])
}

ChannelList2Arr <- function(chnl.lst) {
  stopifnot(is.list(chnl.lst))
  dims <- purrr::map(chnl.lst, dim)
  if (!filesstrings::AllEqual(dims)) {
    stop("The dimensions of the elements of chnl.lst must all be the same.")
  }
  if (length(dims[[1]]) == 2) {
    purrr::reduce(chnl.lst, ~ abind::abind(.x, .y, along = 3))
  } else if (length(dims[[1]]) == 3) {
    purrr::reduce(chnl.lst, ~ abind::abind(.x, .y, along = 4)) %>%
      aperm(c(1, 2, 4, 3))
  } else {
    stop("The elements of chnl.lst must be 2- or 3-dimensional.")
  }
}
