
<!-- index.md is generated from index.Rmd. Please edit that file -->

# `nandb` <img src="man/figures/logo.png" align="right" height=140/>

Calculation of molecular *number and brightness* from fluorescence
microscopy image series. The software was published in a [2016
paper](https://doi.org/10.1093/bioinformatics/btx434). The seminal paper
for the technique is [Digman et
al. 2008](https://doi.org/10.1529/biophysj.107.114645). A
[review](https://doi.org/10.1016/j.ymeth.2017.12.001) of the technique
was published in [2017](https://doi.org/10.1016/j.ymeth.2017.12.001).

If you’re not familiar with the *number and brightness* (N\&B)
technique, then you should familiarise yourself with it by reading the
papers mentioned above before continuing with the `nandb` package. The
`nandb` R package is not intended to introduce people to N\&B, it’s for
people who know about N\&B and want to perform N\&B calculations.

If you’re new to R and you’re here because you want to use `nandb`, be
warned that you will need to learn some basic R first. I recommend
reading the short book “Hands On Programming with R” by Grolemund. This
is available for free at <https://rstudio-education.github.io/hopr/>.
That should be enough but if you want further reading, check out “R for
Data Science” which is available for free at <https://r4ds.had.co.nz/>.

This website gives an introduction to the `nandb` R package, assuming
that the reader has a basic level of N\&B and R knowledge.

## Installation

You can install the release version of `nandb` from
[CRAN](https://CRAN.R-project.org/package=nandb) with:

``` r
install.packages("nandb")
```

You can install the (unstable) development version of `nandb` from
[GitHub](https://github.com/rorynolan/nandb/) with:

``` r
devtools::install_github("rorynolan/nandb")
```

I highly recommend using the release version. The dev version is just
for the ultra-curious and should be thought of as unreliable.

## Using `nandb`

There are two ways to use `nandb`.

1.  Interactively in the R session, playing with the image as a numeric
    array, dealing with one image at a time.
2.  In *batch* mode, having the software read TIFFs, perform the N\&B
    calculations and then write the detrended TIFFs to disk when
    detrending is over. This method permits the user to use R as little
    as possible and is better for those who don’t intend to become bon a
    fide R users.

These are discussed in two articles. These articles deal with
brightness; most people use N\&B to calculate oligomeric state and hence
brightness is the interesting quantity. This package also facilitates
number calculations, which are done in the same way, replacing
“brightness” with “number” in function names. For example, the
“number” equivalent of `brightness_timeseries()` is
`number_timeseries()`. These articles will use the “epsilon” definition
of brightness, but you’re free to use the “B” definition if you prefer
it.

1.  [Brightness calculations on single
    images](https://rorynolan.github.io/nandb/articles/single-images.html)
2.  [Brightness calculations on many images in *batch*
    mode](https://rorynolan.github.io/nandb/articles/batch-mode.html)

Both of these articles mention *brightness timeseries*. These are
explained in the short article [Brightness
timeseries](https://rorynolan.github.io/nandb/articles/brightness-timeseries.html).
N\&B timeseries are a very nice feature of `nandb`, automating a common
and otherwise laborious procedure.
