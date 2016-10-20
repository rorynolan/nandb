#' Make a raster plot of a matrix.
#'
#' Given a matrix \code{mat}, make a raster plot of the matrix whereby in the
#' plot, the pixel at \eqn{x = }\code{i}, \eqn{y = }\code{j} has colour based on
#' the value of \code{mat[i, j]} and the \eqn{x} axis points right and the
#' \eqn{y} axis points down (see "Details").
#'
#' The pixel at \eqn{x = }\code{i}, \eqn{y = }\code{j} has colour based on the
#' value of \code{mat[i, j]} where the \eqn{x} axis points right and the \eqn{y}
#' axis points down. This is in accordance with how
#' \code{\link[EBImage]{EBImage}} and \code{\link{ReadImageData}} (which wraps
#' EBImage's \code{\link[EBImage]{readImage}}). However, when one prints a
#' matrix in a console (or views it in a program such as excel), the value in
#' position \eqn{x = }\code{i}, \eqn{y = }\code{j} is from \code{mat[j, i]}, so
#' if you're confused about a transposed plot, this is why.
#'
#' @param mat The matrix you wish to plot.
#' @param scale A string. The title of the color scale on the right of the plot.
#' @param limits This gives the user the option to set all values outside a
#'   certain range to their nearest value within that range (if \code{clip =
#'   TRUE}) or to \code{NA} (if \code{clip = FALSE}. For example, to set all
#'   values outside the range [1.5, 2.6) to \code{NA}, use \code{limits = c(1.5,
#'   2.6), clip = FALSE}.
#' @param ranges A numeric vector. If you want specific ranges of values to have
#'   the same color, specify these ranges via an increasing numeric vector. For
#'   example, if you want the ranges 0.5-1.2 and 1.2-3, use \code{ranges =
#'   c(0.5, 1.2, 3)}. If \code{ranges} is specified as a number (a numeric
#'   vector of length 1) \code{n}, this is equivalent to setting ranges to be
#'   \code{n} equal-length intervals within the range of the matrix, i.e. it is
#'   equivalent to setting \code{ranges = seq(min(mat), max(mat), length.out = n
#'   + 1)}. At most one of \code{ranges} and \code{limits} should be set. If
#'   ranges is set, the behaviour for values which are not in any of the ranges
#'   are set by the \code{clip} arguments as in the \code{limits} argument.
#' @param colours If you have set \code{ranges}, here you may specify which
#'   colors you wish to colour each range. It must have the same length as the
#'   number of intervals you have specified in \code{ranges}. It is specified as
#'   a character vector, with the colors specified either as the values in
#'   \code{\link{colors}} or as in the value of the \code{\link{rgb}} function.
#'   Note that this allows the use of \code{\link[grDevices]{rainbow}} and
#'   friends. The default uses \code{\link{topo.colors}}.
#' @param cont.colours If ranges is not set, the colours of the pixels are
#'   determined on a continuous colout scale. This scale will have a "low
#'   colour" and a "high colour" (e.g. the default is low values in blue and
#'   high values in red). Again, colours can be specified either as the values
#'   in \code{\link{colors}} or as in the value of the \code{\link{rgb}}
#'   function. For example, to have low values in green and high values in
#'   black, use \code{cont.colours = c("green", "black")}.
#' @param na.colour Which colour should the \code{NA} pixels be? Default is
#'   black.
#' @param clip If either \code{limits} or \code{ranges} are set (one should
#'   never set both), there may be values that fall outside the specified
#'   limits/ranges. If \code{clip = TRUE}, values outside these limits/ranges
#'   are set to their nearest values within them, but if \code{clip = FALSE},
#'   these values are set to NA. Note that setting \code{clip = TRUE} is
#'   equivalent to setting both \code{clip.low} and \code{clip.high} to
#'   \code{TRUE}.
#' @param clip.low Setting this to \code{TRUE} (and leaving \code{clip = FALSE},
#'   \code{clip.high = FALSE}) will set all values falling below the specified
#'   limits/ranges to their nearest value within them, but all values falling
#'   above those limits/ranges will be set to \code{NA}.
#' @param clip.high Setting this to \code{TRUE} (and leaving \code{clip =
#'   FALSE}, \code{clip.low = FALSE}) will set all values falling above the
#'   specified limits/ranges to their nearest value within them, but all values
#'   falling below those limits/ranges will be set to \code{NA}.
#'
#' @return In the graphics console, a raster plot (via
#'   \code{\link[ggplot2]{geom_raster}}) will appear with the matrix values
#'   represented as pixel colours, with a named scale bar.
#' @export
MatrixRasterPlot <- function(mat, scale.name = "scale",
                             limits = NULL, ranges = NULL,
                             colours = NULL, cont.colours = c("blue", "red"),
                             na.colour = "black", clip = FALSE,
                             clip.low = FALSE, clip.high = FALSE) {
  plain.theme <- ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                axis.text = ggplot2::element_blank(),
                                axis.title = ggplot2::element_blank(),
                                panel.background = ggplot2::element_rect(fill = "white"))
  df <- reshape2::melt(mat) %>%
    dplyr::transmute(x = Var1, y = 1 + max(Var2) - Var2, value = value)
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
    df <- df %>% dplyr::mutate(colour = colours.ranges)
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = colour)) +
      ggplot2::scale_fill_manual(scale.name, values = topocols) +
      ggplot2::geom_raster() +
      plain.theme
  } else {
    if (is.null(limits)) limits <- range(df$value)
    if (clip) {
      clip.low <- TRUE
      clip.high <- TRUE
    }
    if (clip.low) df$value[df$value < limits[1]] <- limits[1]
    if (clip.high) df$value[df$value > limits[2]] <- limits[2]
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::scale_fill_gradient(scale.name, limits = limits,
                                   low = cont.colours[1], high = cont.colours[2],
                                   na.value = na.colour) +
      ggplot2::geom_raster() +
      plain.theme
  }
}

