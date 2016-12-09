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
#'   2.6), clip = FALSE}. The colour range will cover all values within these
#'   specified limits.
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
#'   number of intervals you have specified in \code{ranges}.  If you have not
#'   specified \code{ranges}, here you may specify the colours (to be passed to
#'   \code{\link[ggplot2]{scale_fill_gradientn}}) to create the continuous
#'   colour band. It is specified as a character vector, with the colors
#'   specified either as the values in \code{\link{colors}} or as in the value
#'   of the \code{\link{rgb}} function. Note that this allows the use of
#'   \code{\link[grDevices]{rainbow}} and friends. The default uses
#'   \code{\link[viridis]{viridis}}.
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
#' @param log.trans Do you want to log-transform the colour scaling?
#' @param breaks Where do you want tick marks to appear on the legend colour scale?
#'
#' @return In the graphics console, a raster plot (via
#'   \code{\link[ggplot2]{geom_raster}}) will appear with the matrix values
#'   represented as pixel colours, with a named scale bar.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file("extdata",
#' "low_oligomers.tif",
#' package = "nandb"))
#' display(normalize(img[, , 1]), method = "raster")
#' brightness <- Brightness(img, tau = 10, mst = "Huang")
#' MatrixRasterPlot(brightness, scale.name = "brightness")
#' MatrixRasterPlot(brightness, scale.name = "brightness", log.trans = TRUE)
#' MatrixRasterPlot(brightness, scale.name = "brightness", log.trans = TRUE,
#' breaks = 1:3)
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#' ranges = seq(0.5, 3, length.out = 6), range.names = paste0(1:5, "mer"))
#' @export
MatrixRasterPlot <- function(mat, scale.name = "scale",
                             limits = NULL, ranges = NULL,
                             range.names = NULL, colours = NULL,
                             na.colour = "black", clip = FALSE,
                             clip.low = FALSE, clip.high = FALSE,
                             log.trans = FALSE, breaks = NULL) {
  plain.theme <- ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                axis.text = ggplot2::element_blank(),
                                axis.title = ggplot2::element_blank(),
                                panel.background = ggplot2::element_rect(fill = "white"),
                                legend.key.height = unit(1, "cm"))
  rownames(mat) <- NULL
  colnames(mat) <- NULL
  df <- reshape2::melt(mat) %>%
    dplyr::transmute(x = Var1, y = 1 + max(Var2) - Var2, value = value)
  if (!is.null(ranges)) {
    nr <- length(ranges - 1)
    if (is.null(colours)) {
      if (log.trans) {
        min.log.sep <- log(nr) - log(nr - 1)
        nums <- round((1 + log(seq_along(ranges))) / min.log.sep)
        colours <- viridis::viridis(max(nums))[nums]
      } else {
        colours <- viridis::viridis(nr)
      }
    } else if (length(colours) != nr) {
      stop("The number of colours must match the number of ranges")
    }
    ranges <- ranges %>% {cbind(.[-length(.)], .[-1])}  # adjacent pairs
    ranges.typeset <- apply(ranges, 1,
                            function(x) paste(round(x, 2), collapse = "-"))
    colours.ranges <- factor(df$value %>% sapply(WhichInterval, ranges),
                             levels = seq_len(nrow(ranges)))
    df <- dplyr::mutate(df, colour = colours.ranges)
    if (!is.null(range.names)) {
      if (length(range.names) != nrow(ranges)) {
        stop("The number of range.names must be equal to the number of ranges (specified directly via ranges or indirectly via colours.")
      }
    } else {
      range.names <- levels(colours.ranges)
    }
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = colour)) +
      ggplot2::scale_fill_manual(scale.name,
                                 values = colours %T>% {
                                   names(.) = as.character(seq_along(colours))
                                 },
                                 na.value = na.colour,
                                 labels = range.names %T>% {
                                   names(.) = as.character(seq_along(colours))
                                 }) +
      ggplot2::geom_raster() +
      plain.theme + coord_fixed()
  } else {
    if (is.null(limits)) limits <- range(df$value)
    if (is.null(colours)) colours <- viridis::viridis(9)
    if (clip) {
      clip.low <- TRUE
      clip.high <- TRUE
    }
    if (clip.low) df$value[df$value < limits[1]] <- limits[1]
    if (clip.high) df$value[df$value > limits[2]] <- limits[2]
    if (is.null(breaks)) {
      if (log.trans) {
        if (min(limits) <= 0) {
          stop("The data or limits range are not all positive-valued therefore a log-transformation is not possible.")
        } else {
          # limits <- log(limits)
        }
        breaks <- function(x) {
          untruncated <- exp(seq(log(min(x)), log(max(x)), length.out = 5))
          min.sep <- min(diff(untruncated))
          roundn <- log10(min.sep) %>% {  # get a sensible amount of digits for breaks
            ifelse(. >= 1, 0, -floor(.))
          }
          round(untruncated, roundn)
        }
      } else {
        breaks <- waiver()
      }
    }
    ggplot2::ggplot(df, ggplot2::aes(x, y, fill = value)) +
      ggplot2::scale_fill_gradientn(scale.name, limits = limits,
                                   colours = colours,
                                   na.value = na.colour,
                                   trans = ifelse(log.trans,
                                                  "log", "identity"),
                                   breaks = breaks) +
      ggplot2::geom_raster() + plain.theme + coord_fixed()
  }
}

#' A brightness image with a different colour for each kmer.
#'
#' Make a colour image based on a brightness image where each kmer has its own
#' colour. This requires you te specify the brightness of a monomer (which
#' should be greater than 1).
#'
#' @param brightness.mat The brightness matrix.
#' @param monomer.brightness The (median) brightness of a monomer.
#' @param log.trans Do you want to log-transform the colour scaling?
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file("extdata",
#' "FKBP-mClover_before_0.5nM_AP1510.tif",
#' package = "nandb"))
#' display(normalize(img[, , 1]), method = "raster")
#' brightness <- Brightness(img, tau = 10, mst = "Huang")
#' KmerPlot(brightness, 1.16)
#'
#' @export
KmerPlot <- function(brightness.mat, monomer.brightness, log.trans = FALSE) {
  stopifnot(is.matrix(brightness.mat))
  stopifnot(monomer.brightness > 1)
  max.b <- max(brightness.mat, na.rm = TRUE)
  if (max.b > monomer.brightness) {
    ranges <- c(0, seq(1 + 0.5 * (monomer.brightness - 1),
                       max.b, monomer.brightness - 1),
                max.b) %>% unique  # the unique avoids the unlikely possibility of repeating the max at the end
    lrm2 <- length(ranges) - 2  # LengthofRangesMinus2
    range.names <- c("Immobile", paste0(seq_len(lrm2), "mers"))
    if (log.trans) {
      min.log.sep <- log(lrm2) - log(lrm2 - 1)
      nums <- round((1 + log(seq_len(lrm2))) / min.log.sep)
      colours <- c("slategray4", viridis::viridis(max(nums))[nums])
    } else {
      colours <- c("slategray4", viridis::viridis(lrm2))
    }
  } else {
    ranges <- c(0, max.b)
    range.names <- "Immobile"
    colours = "slategray4"
  }
  MatrixRasterPlot(brightness.mat, scale.name = "Brightness", ranges = ranges,
                   range.names = range.names, colours = colours,
                   log.trans = log.trans)
}

#' Make brightness plots (images) for an entire folder.
#'
#' This requires that the folder already have the brightnesses saved as csv
#' files. Output images are saved as pdf.
#'
#' @param folder.path The path (relative or absolute) to the folder you wish to
#'   process.
#' @param patt The pattern (in regular expression) of the csv files in the
#'   folder that are brightness images. The default matches any .csv files which
#'   contains "Brightness" or "brightness" in its name.
#' @param verbose Get a real time report on progress?
#' @param ... Parameters to pass to \code{\link{MatrixrasterPlot}} (these
#'   \emph{must} be named arguments).
#'
BrightnessPlotFolder <- function(folder.path = ".",
                                 patt = "[Bb]rightness.*\\.csv$",
                                 verbose = TRUE,
                                 ...) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  file.names <- list.files(pattern = patt)
  brightness.matrices <- lapply(file.names, ReadImageTxt)
  dots <- list(...)
  if (is.null(dots$limits)) {
    dots$limits <- range(unlist(brightness.matrices), na.rm = TRUE)
  }
  if (is.null(dots$scale.name)) dots$scale.name <- "brightness"
  for (i in seq_along(file.names)) {
    bld <- filesstrings::BeforeLastDot(file.names[i])
    if (verbose) {
      paste0("Now processing ", bld, ".") %>% message
    }
    pdf.file.name <- bld %>% filesstrings::MakeExtName("pdf")
    pdf(pdf.file.name)
    c(list(mat = brightness.matrices[[i]]), dots) %>%
    do.call(MatrixRasterPlot, .) %>% print
    dev.off()
  }
}
