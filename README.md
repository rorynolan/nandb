nandb
================

An R package for performing number and brightness analysis as in Digman et al. 2008.

[![Travis-CI Build Status](https://travis-ci.org/rorynolan/nandb.svg?branch=master)](https://travis-ci.org/rorynolan/nandb) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rorynolan/nandb?branch=master&svg=true)](https://ci.appveyor.com/project/rorynolan/nandb) [![codecov](https://codecov.io/gh/rorynolan/nandb/branch/master/graph/badge.svg)](https://codecov.io/gh/rorynolan/nandb) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/nandb)](https://cran.r-project.org/package=nandb) ![RStudio CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/nandb) ![RStudio CRAN monthly downloads](http://cranlogs.r-pkg.org/badges/nandb) [![Rdocumentation](http://www.rdocumentation.org/badges/version/nandb)](http://www.rdocumentation.org/packages/nandb) ![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)

Installation
------------

### Linux

1.  On **Ubuntu** (similarly for other debian **linux**), you need to do:

        sudo apt-get update
        sudo apt-get install libssl-dev libtiff5-dev libfftw3-dev 
        sudo apt-get install libcurl4-openssl-dev libxml2-dev 

2.  In R, run

``` r
install.packages("nandb")
```

### Windows

1.  Go to <https://cran.r-project.org/bin/windows/Rtools/> and install the latest version of Rtools.

2.  In R, run

``` r
install.packages("nandb")
```

### Mac

In R, run

``` r
install.packages("nandb")
```

Installing the development version
----------------------------------

The release version is recommended (and installed with `install.packages("nandb")` as above), but to install the development version, in R enter

``` r
devtools::install_github("rorynolan/nandb")
```

Use
---

For the lowdown on how to use this package, you should read the package [vignette](https://cran.rstudio.com/web/packages/nandb/vignettes/nandb.html). Once you've installed the package, you can browse them using `vignette(package = "nandb")`.

Contribution
------------

Contributions to this package are welcome. The preferred method of contribution is through a github pull request. Feel free to contact me by creating an issue. Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
