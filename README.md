nandb
================

An R package for performing number and brightness analysis as in Digman et al. 2008.

[![Travis-CI Build Status](https://travis-ci.org/rorynolan/nandb.svg?branch=master)](https://travis-ci.org/rorynolan/nandb) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rorynolan/nandb?branch=master&svg=true)](https://ci.appveyor.com/project/rorynolan/nandb) [![codecov](https://codecov.io/gh/rorynolan/nandb/branch/master/graph/badge.svg)](https://codecov.io/gh/rorynolan/nandb) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/nandb)](https://cran.r-project.org/package=nandb) ![RStudio CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/nandb) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

Installation
------------

### Platform Dependencies

#### Java Stuff

##### Everyone except Mac OS X

Install the latest [JDK](http://www.oracle.com/technetwork/java/javase/downloads/index.html).

##### Mac OS X

Mac OS X comes with a legacy Apple Java 6. Update your Java installation to a newer version provided by Oracle.

1.  Install [Oracle JDK](http://www.oracle.com/technetwork/java/javase/downloads/index.html).

2.  Update R Java configuration by executing from the command line (you might have to run it as a super user by prepending `sudo` depending on your installation).

        R CMD javareconf

3.  Re-install *rJava* from sources in order to properly link to the non-system Java installation.

    ``` r
    install.packages("rJava", type="source")
    ```

Then, to make everything work with Rstudio, run

    sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib

(taken from <http://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite>).

You can verify your configuration by running the following commands. This should return the Java version string corresponding to the one downloaded and installed in step 1.

``` r
library(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
## [1] "1.8.0_112-b16" 
```

Thanks to @aoles for these instructions. Check out his great work on image analysis in R at <https://github.com/aoles>.

##### Other stuff

This is slightly painful but you'll only have to do it once. First of all, if you're on **Ubuntu** (similarly for other debian **linux**), you need to do:

    sudo apt-get update
    sudo apt-get install libssl-dev libtiff5-dev libfftw3-dev 
    sudo apt-get install libcurl4-openssl-dev libxml2-dev 
    sudo apt-get install default-jre default-jdk libboost-all-dev

On **mac**, you need to go to <https://support.apple.com/kb/dl1572?locale=en_US> and download and install Java 6 (don't ask me why this is necessary). Then you need to install `Xcode` from the app store. Finally, you need to open a terminal and type `brew install boost`. If `brew install boost` doesn't work then google "install boost c++ library mac"; you need to have this done in order to use `nandb`.

On **Windows**, you need to go to <https://cran.r-project.org/bin/windows/Rtools/> and install the latest version of Rtools.

### All Platforms

Then, everyone, open R and run:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("EBImage", "BiocParallel"))
install.packages("nandb")
```

If you get a message saying `Update all/some/none? [a/s/n]:`, type `a`. If you get the chance to install some packages from source, it's safest to decline.

#### Problems

If you run into problems during your installation, it's most likely that your installation of `rJava` didn't work. Try running `install.packages("rJava")` and try to work through the errors there. If that still doesn't work, try googling "install rJava" for your operating system e.g. "install rJava Ubuntu". The second most likely culprit is `EBImage` so similarly have a google of "install EBImage". If you get these to work, then try the installation instructions again. If you still can't get it to work, feel free to contact me via the issues page associated with this repo.

### Updates

To update the package, you just need to do the same thing:

``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("EBImage", "BiocParallel"))
install.packages("nandb")
```

Use
---

For the lowdown on how to use this package, you should read the package vignettes. Once you've installed the package, you can browse them using `vignette(package = "nandb")`.
