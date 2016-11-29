nandb
================

An R package for performing number and brightness analysis as in Digman et al. 2008.

Installation
------------

### Platform Dependencies

This is slightly painful but you'll only have to do it once. First of all, if you're on **Ubuntu** (similarly for other debian **linux**), you need to do:

    sudo apt-get update
    sudo apt-get install libssl-dev libtiff5-dev libfftw3-dev 
    sudo apt-get install libcurl4-openssl-dev libxml2-dev 
    sudo apt-get install default-jre default-jdk

On **mac**, you need to go to <https://support.apple.com/kb/dl1572?locale=en_US> and download and install Java 6 (don't ask me why this is necessary).

On **Windows**, you need to go to <https://cran.r-project.org/bin/windows/Rtools/> and install the latest version of Rtools.

### All Platforms

Then, everyone, open R and run:

``` r
install.packages(c("devtools", "knitr", "rmarkdown", "XML", "xml2"))
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

If you get a message saying `Update all/some/none? [a/s/n]:`, type `a`. Then run

``` r
options(unzip = "internal")
library(devtools)
install_github("rorynolan/filesstrings")
install_github("rorynolan/autothresholdr")
install_github("rorynolan/nandb",
               auth_token = "50627d2ca9badc802ed993cb2b9a914034ab2246",
               build_vignettes = TRUE)
```

Done.

#### Problems

If you run into problems during your installation, it's most likely that your installation of `rJava` didn't work. Try running `install.packages("rJava")` and try to work through the errors there. If that still doesn't work, try googling "install rJava" for your operating system e.g. "install rJava Ubuntu". The second most likely culprit is `EBImage` so similarly have a google of "install EBImage". If you get these to work, then try the installation instructions again. If you still can't get it to work, feel free to contact us via the issues page associated with this repo.

### Updates

To update the package, you just need to do the same thing:

``` r
options(unzip = "internal")
library(devtools)
install_github("rorynolan/filesstrings")
install_github("rorynolan/autothresholdr")
install_github("rorynolan/nandb",
               auth_token = "50627d2ca9badc802ed993cb2b9a914034ab2246",
               build_vignettes = TRUE)
```

To check if you need an update, check if the package has been updated since you installed it. To check your current version, use `packageVersion("nandb")`. To check if there's a newer version, go to the github page \[github.com/rorynolan/filesstrings\] (you're probably there right now) and look for the version in the description file.

If you don't want to bother checking and you just want to make sure you have the latest version, just run those two lines of code.

Use
---

For the lowdown on how to use this package, you should read the package vignette (using `vignette("nandb")`); however, here's a quick example. First load the libraries:

``` r
library("nandb")
library("magrittr")
```

Then load the .tif file included with the package:

``` r
tif_file_location <- system.file("extdata", 
                                 "FKBP-mClover_before_0.5nM_AP1510.tif", 
                                 package = "nandb")
img <- ReadImageData(tif_file_location)
```

This `img` is a 3-dimensional array of image slices, where each element corresponds to a pixel with a value of integer counts. Looking at the dimensions of this object:

``` r
dim(img)
#> [1] 256 256  50
```

we can see that it is 50 slices each of which is a 256x256 image. We can view the first slice:

``` r
EBImage::display(EBImage::normalize(img[, , 1]), method = "raster")
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Now we can calculate the brightness image based on this image series, with an exponential filtering detrend with a time constant of (say) 10 frames via

``` r
img_brightness <- Brightness(img, tau = 10, mst = "Huang")
```

Now we can see the (median filtered) version of the output.

``` r
img_brightness %>% MedianFilterB %>% 
  MatrixRasterPlot(scale.name = "brightness", limits = c(0.9, 1.5))
```

![](README_files/figure-markdown_github/Brightness%20plot-1.png)
