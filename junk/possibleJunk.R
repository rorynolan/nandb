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
