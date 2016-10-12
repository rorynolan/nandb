nandb
================

An R package for performing number and brightness analysis as in Digman et al. 2008.

Installation
------------

This is slightly painful but you'll only have to do it once. Open R and run:

``` r
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

If you get a message saying `Update all/some/none? [a/s/n]:`, type `a`. Then run

``` r
devtools::install_github("rorynolan/filesstrings")
devtools::install_github("rorynolan/nandb")
```

Done.

Updates
-------

To update the package, you just need

``` r
devtools::install_github("rorynolan/filesstrings")
devtools::install_github("rorynolan/nandb")
```

To check if you need an update, check if the package has been updated since you installed it. To check your current version, use `packageVersion("nandb")`. To check if there's a newer version, go to the github page \[github.com/rorynolan/filesstrings\] (you're probably there right now) and look for the version in the description file.

If you don't want to bother checking and you just want to make sure you have the latest version, just run those two lines of code.
