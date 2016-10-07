MatrixPlot <- function(mat, scale.name = "scale", ranges = NULL,
                       colours = NULL, limits = NULL, low = "blue", high = "red",
                       clip = FALSE, clip.low = FALSE, clip.high = FALSE) {
  plain.theme <- theme(axis.ticks = element_blank(), axis.text = element_blank(),
                       axis.title = element_blank(),
                       panel.background = element_rect(fill = "white"))
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
      scale_fill_manual(scale.name, values = topocols) +
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
      scale_fill_gradient(scale.name, limits = limits, low = low, high = high, na.value = "black") +
      geom_raster() +
      plain.theme
  }
}
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
