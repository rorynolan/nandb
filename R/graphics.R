#' Make a raster plot of a matrix.
#'
#' Given a matrix `mat`, make a raster plot of the matrix whereby in the
#' plot, the pixel at \eqn{x = }`i`, \eqn{y = }`j` has colour based on
#' the value of `mat[i, j]` and the \eqn{x} axis points right and the
#' \eqn{y} axis points down (see 'Details').
#'
#' The pixel at \eqn{x = }`i`, \eqn{y = }`j` has colour based on the
#' value of `mat[i, j]` where the \eqn{x} axis points right and the \eqn{y}
#' axis points down. This is in accordance with how
#' [EBImage::EBImage()] and [ReadImageData()] (which wraps
#' EBImage's [EBImage::readImage()]). However, when one prints a
#' matrix in a console (or views it in a program such as excel), the value in
#' position \eqn{x = }`i`, \eqn{y = }`j` is from `mat[j, i]`, so
#' if you're confused about a transposed plot, this is why.
#'
#' @param mat The matrix you wish to plot.
#' @param scale.name A string. The title of the color scale on the right of the
#'   plot.
#' @param limits This gives the user the option to set all values outside a
#'   certain range to their nearest value within that range (if \code{clip =
#'   TRUE}) or to `NA` (if `clip = FALSE`. For example, to set all
#'   values outside the range [1.5, 2.6) to `NA`, use \code{limits = c(1.5,
#'   2.6), clip = FALSE}. The colour range will cover all values within these
#'   specified limits.
#' @param ranges A numeric vector. If you want specific ranges of values to have
#'   the same color, specify these ranges via an increasing numeric vector. For
#'   example, if you want the ranges 0.5-1.2 and 1.2-3, use \code{ranges =
#'   c(0.5, 1.2, 3)}. If `ranges` is specified as a number (a numeric
#'   vector of length 1) `n`, this is equivalent to setting ranges to be
#'   `n` equal-length intervals within the range of the matrix, i.e. it is
#'   equivalent to setting \code{ranges = seq(min(mat), max(mat), length.out = n
#'   + 1)}. At most one of `ranges` and `limits` should be set. If
#'   ranges is set, the behaviour for values which are not in any of the ranges
#'   are set by the `clip` arguments as in the `limits` argument.
#' @param range.names A character vector. If your colour scale is discrete, here
#'   you can set the names which will label each range in the legend.
#' @param colours If you have set `ranges`, here you may specify which
#'   colors you wish to colour each range. It must have the same length as the
#'   number of intervals you have specified in `ranges`.  If you have not
#'   specified `ranges`, here you may specify the colours (to be passed to
#'   [ggplot2::scale_fill_gradientn()]) to create the continuous
#'   colour band. It is specified as a character vector, with the colors
#'   specified either as the values in [colors()] or as in the value
#'   of the [rgb()] function. Note that this allows the use of
#'   [grDevices::rainbow()] and friends. The default uses
#'   [viridis::viridis()].
#' @param na.colour Which colour should the `NA` pixels be? Default is
#'   black.
#' @param clip If either `limits` or `ranges` are set (one should
#'   never set both), there may be values that fall outside the specified
#'   limits/ranges. If `clip = TRUE`, values outside these limits/ranges
#'   are set to their nearest values within them, but if `clip = FALSE`,
#'   these values are set to NA. Note that setting `clip = TRUE` is
#'   equivalent to setting both `clip.low` and `clip.high` to
#'   `TRUE`.
#' @param clip.low Setting this to `TRUE` (and leaving `clip = FALSE`,
#'   `clip.high = FALSE`) will set all values falling below the specified
#'   limits/ranges to their nearest value within them, but all values falling
#'   above those limits/ranges will be set to `NA`.
#' @param clip.high Setting this to `TRUE` (and leaving \code{clip =
#'   FALSE}, `clip.low = FALSE`) will set all values falling above the
#'   specified limits/ranges to their nearest value within them, but all values
#'   falling below those limits/ranges will be set to `NA`.
#' @param log.trans Do you want to log-transform the colour scaling?
#' @param breaks Where do you want tick marks to appear on the legend colour
#'   scale?
#' @param include.breaks If you don't want to specify all the breaks, but you
#'   want some specific ones to be included on the legend colour scale, specify
#'   those specific ones here.
#'
#' @return In the graphics console, a raster plot (via
#'   [ggplot2::geom_raster()]) will appear with the matrix values
#'   represented as pixel colours, with a named scale bar.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' brightness <- Brightness(img, tau = NA, mst = "Huang")
#' MatrixRasterPlot(brightness, scale.name = 'brightness')
#' MatrixRasterPlot(brightness, scale.name = 'brightness', log.trans = TRUE)
#' MatrixRasterPlot(brightness, scale.name = 'brightness', log.trans = TRUE,
#'                  include.breaks = 1.35)
#' MatrixRasterPlot(brightness, scale.name = 'brightness', log.trans = TRUE,
#'                  breaks = 1:3)
#' MatrixRasterPlot(brightness, scale.name = 'brightness',
#'   ranges = seq(0.5, 3, length.out = 6), range.names = paste0(1:5, 'mer'))
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#'                  ranges = seq(0.5, 3, length.out = 6),
#'                  range.names = paste0(1:5, "mer"), log.trans = TRUE)
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#'   include.breaks = 1.25, range.names = NULL, log.trans = FALSE)
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#'                  include.breaks = 1.25, log.trans = TRUE)
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#'                  limits = c(1, 1.25), clip = TRUE)
#' MatrixRasterPlot(brightness, scale.name = "brightness",
#'                  include.breaks = 1.25)
#'
#' @export
MatrixRasterPlot <- function(mat, scale.name = "scale", limits = NULL,
  ranges = NULL, range.names = NULL, colours = NULL, na.colour = "black",
  clip = FALSE, clip.low = FALSE, clip.high = FALSE, log.trans = FALSE,
  breaks = NULL, include.breaks = NULL) {
  stopifnot(is.matrix(mat))
  plain.theme <- theme(axis.ticks = element_blank(),
    axis.text = element_blank(), axis.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.key.height = unit(1, "cm"))
  rownames(mat) <- NULL
  colnames(mat) <- NULL
  include.breaks <- unique(include.breaks)
  df <- reshape2::melt(mat) %>% dplyr::transmute(x = Var1,
    y = 1 + max(Var2) - Var2, value = value)
  if (!is.null(ranges)) {
    nr <- length(ranges) - 1
    if (is.null(colours)) {
      if (log.trans) {
        min.log.sep <- log(nr) - log(nr - 1)
        nums <- round((1 + log(seq_along(ranges))) / min.log.sep) %>%
          {
          .[-length(.)]
          }
        colours <- viridis::viridis(max(nums))[nums]
      } else {
        colours <- viridis::viridis(nr)
      }
    } else if (length(colours) != nr) {
      stop("The number of colours must match the number of ranges")
    }
    ranges <- ranges %>% {
      cbind(.[-length(.)], .[-1])
    }  # adjacent pairs
    ranges.typeset <- apply(ranges, 1, function(x) paste(round(x,
      2), collapse = "-"))
    colours.ranges <- factor(df$value %>% vapply(WhichInterval, integer(1),
      ranges), levels = seq_len(nrow(ranges)))
    df <- dplyr::mutate(df, colour = colours.ranges)
    if (!is.null(range.names)) {
      if (length(range.names) != nrow(ranges)) {
        stop("The number of range.names must be equal to the number of ranges ",
             "(specified directly via ranges or indirectly via colours.")
      }
    } else {
      range.names <- levels(colours.ranges)
    }
    ggplot(df, aes(x, y, fill = colour)) + scale_fill_manual(scale.name,
      values = magrittr::set_names(colours, seq_along(colours)),
      na.value = na.colour,
      labels = magrittr::set_names(range.names, seq_along(colours))) +
      geom_raster() + plain.theme + coord_fixed()
  } else {
    if (is.null(limits)) {
      if (is.null(include.breaks)) {
        also <- numeric(0)
      } else {
        also <- include.breaks
      }
      limits <- range(c(df$value, also), na.rm = TRUE)
    } else {
      if (!is.null(include.breaks))
        limits <- range(c(limits, include.breaks))
    }
    if (is.null(colours)) colours <- viridis::viridis(99)
    if (clip) {
      clip.low <- TRUE
      clip.high <- TRUE
    }
    if (clip.low) df$value[df$value < limits[1]] <- limits[1]
    if (clip.high) df$value[df$value > limits[2]] <- limits[2]
    if (is.null(breaks)) {
      if (log.trans) {
        if (min(limits) <= 0) {
          stop("The data or limits range are not all positive-valued, ",
          "therefore a log-transformation is not possible.")
        }
        if (is.null(include.breaks)) {
          breaks <- ScaleBreaksFun(include.breaks = include.breaks,
                                   log = log.trans)
        } else {
          breaks <- ScaleBreaksFun(include.breaks = include.breaks,
                                   log = log.trans)
        }
      } else {
        if (is.null(include.breaks)) {
          breaks <- waiver()
        } else {
          breaks <- ScaleBreaksFun(include.breaks = include.breaks,
                                   log = log.trans)
        }
      }
    }
    ggplot(df, aes(x, y, fill = value)) + scale_fill_gradientn(scale.name,
      limits = limits, colours = colours, na.value = na.colour,
      trans = ifelse(log.trans, "log", "identity"), breaks = breaks) +
      geom_raster() + plain.theme + coord_fixed()
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
#' This is a `ggplot2` object and can be manipulated thus.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' brightness <- Brightness(img, tau = NA, mst = "Huang")
#' KmerPlot(brightness, 1.02)
#' KmerPlot(brightness, 1.12)
#' KmerPlot(brightness, 100)
#' KmerPlot(brightness, 1.02, log.trans = TRUE)
#' KmerPlot(MedianFilterB(brightness), 1.02)
#'
#' @export
KmerPlot <- function(brightness.mat, monomer.brightness, log.trans = FALSE) {
  stopifnot(is.matrix(brightness.mat))
  stopifnot(monomer.brightness > 1)
  max.b <- max(brightness.mat, na.rm = TRUE)
  if (max.b > monomer.brightness) {
    ranges <- c(0, seq(1 + 0.5 * (monomer.brightness - 1),
      max.b, monomer.brightness - 1), max.b) %>% unique %>% {
        # the unique avoids the unlikely possibility of repeating the
        # max at the end
        if (length(.) <= 13) {
          .
        } else {
          CollapseRanges(., 13, prefer.low = TRUE)
        }
      }
    lrm2 <- length(ranges) - 2  # LengthofRangesMinus2
    range.names <- c("Immobile", paste0(seq_len(lrm2), "mers"))
    range.names[length(range.names)] <- paste0(">=",
                                               range.names[length(range.names)])
    if (log.trans) {
      min.log.sep <- log(lrm2) - log(lrm2 - 1)
      nums <- round((1 + log(seq_len(lrm2)))/min.log.sep)
      colours <- c("slategray4", viridis::viridis(max(nums))[nums])
    } else {
      colours <- c("slategray4", viridis::viridis(lrm2))
    }
  } else {
    ranges <- c(0, max.b)
    range.names <- "Immobile"
    colours <- "slategray4"
  }
  MatrixRasterPlot(brightness.mat, scale.name = "Brightness",
    ranges = ranges, range.names = range.names, colours = colours,
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
#'   contains 'Brightness' or 'brightness' in its name.
#' @param verbose Get a real time report on progress?
#' @param ... Parameters to pass to [MatrixRasterPlot()] (these
#'   \emph{must} be named arguments).
#'
#' @return This is a `ggplot2` object and can be manipulated thus.
#'
#' @examples
#' setwd(tempdir())
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' WriteIntImage(img, '50.tif')
#' WriteIntImage(img, '50again.tif')
#' BrightnessTxtFolder(tau = NA, mst = "Huang", mcc = 2)
#' BrightnessPlotFolder()
#' list.files()
#' file.remove(list.files())  # cleanup
#'
#' @export
BrightnessPlotFolder <- function(folder.path = ".",
  patt = ".*[Bb]rightness.*\\.csv$", verbose = TRUE, ...) {
  init.dir <- getwd()
  on.exit(setwd(init.dir))
  setwd(folder.path)
  file.names <- list.files(pattern = patt)
  brightness.matrices <- lapply(file.names, ReadImageTxt)
  dots <- list(...)
  if (is.null(dots$limits)) {
    dots$limits <- range(unlist(brightness.matrices), na.rm = TRUE)
  }
  if (is.null(dots$scale.name))
    dots$scale.name <- "brightness"
  for (i in seq_along(file.names)) {
    bld <- filesstrings::before_last_dot(file.names[i])
    if (verbose) {
      paste0("Now processing ", bld, ".") %>% message
    }
    pdf.file.name <- filesstrings::give_ext(bld, "pdf")
    grDevices::pdf(pdf.file.name)
    c(list(mat = brightness.matrices[[i]]), dots) %>% do.call(MatrixRasterPlot,
      .) %>% print
    grDevices::dev.off()
  }
}

#' Plot the values in two arrays against each other.
#'
#' Plot the values of two arrays of identical dimension against each other using
#' a hexagonal heatmap.
#'
#' @param arr.x,arr.y The two arrays. The `arrx` values will be along the
#'   \eqn{x} axis and the `arry` values along the \eqn{y} axis.
#' @param bins Passed to [ggplot2::geom_hex()].
#' @param log.trans Do you want to log-transform the colour scaling?
#' @param colours Here you may specify the colours (to be passed to
#'   [ggplot2::scale_fill_gradientn()]) to create the continuous colour band. It
#'   is specified as a character vector, with the colors specified either as the
#'   values in [colors()] or as in the value of the [rgb()] function. Note that
#'   this allows the use of [grDevices::rainbow()] and friends. The default uses
#'   [viridis::viridis()].
#' @param limits A numeric vector of length two providing limits of the scale.
#' @param breaks Where do you want tick marks to appear on the legend colour
#'   scale?
#' @param include.breaks If you don't want to specify all the breaks, but you
#'   want some specific ones to be included on the legend colour scale, specify
#'   those specific ones here.
#'
#' @return This is a `ggplot2` object and can be manipulated thus.
#'
#' @examples
#' library(EBImage)
#' img <- ReadImageData(system.file('extdata', '50.tif', package = 'nandb'))
#' display(normalize(img[, , 1]), method = 'raster')
#' brightness <- Brightness(img, tau = NA, mst = "Huang")
#' mean.intensity <- MeanPillars(img)
#' ArrArrHexPlot(mean.intensity, brightness) +
#' ggplot2::labs(x = 'intensity', y = 'brightness')
#'
#' @export
ArrArrHexPlot <- function(arr.x, arr.y, bins = 60, log.trans = FALSE,
  colours = NULL, limits = NULL, breaks = NULL, include.breaks = NULL) {
  stopifnot(all.equal(dim(arr.x), dim(arr.y)))
  df <- data.frame(x = as.vector(arr.x), y = as.vector(arr.y))
  if (is.null(colours))
    colours <- viridis::viridis(9)
  if (is.null(limits))
    limits <- NULL
  if (is.null(breaks)) breaks <- ScaleBreaksFun(include.breaks, log.trans)
  ggplot(df, aes(x, y)) + geom_hex(bins = bins) +
    scale_fill_gradientn(colours = colours, breaks = breaks,
    limits = limits, trans = ifelse(log.trans, "log", "identity"))
}
