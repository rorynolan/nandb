## 1.0.0

#### NEW FEATURES
* The package is now peer-reviewed and published in the journal *Bioinformatics*. See `citation("nandb")`.
* The package style is now in accordance with the tidyverse style guide.
* `brightness()` and `number()` now include options to set `offset`, S-factor, `readout_noise` and `gamma` correction terms.
* `brightness()` and `number()` now enable calculation of both definitions ("B" and "epsilon"; "N" and "n") of brightness and number.
* Detrending and tiff I/O are outsourced to the `detrendr` package. This makes detrending more accurate and much faster.
* The package now has its own S3 class system.
* Graphics are streamlined to one `display()` function.


#### BUG FIXES 
* Kmer calculations are no longer possible. The way in which they were done was over-simple.


### 0.2.1

#### BUG FIXES
* Compatible with `filesstrings` 1.1.0.


## nandb 0.2.0

* The first version that I consider CRAN-worthy.



