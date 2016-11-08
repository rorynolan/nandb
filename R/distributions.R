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
  dens <- density(samp)
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
