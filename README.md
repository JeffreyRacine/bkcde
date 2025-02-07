# bkcde

This is the R package `bkcde' (Boundary Corrected Polynomial Adaptive Univariate Conditional Density Estimation) written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca)

## Installation

You can install the more recent version from GitHub by downloading the [zip ball](https://github.com/JeffreyRacine/R-Package-bkcde/zipball/main) or [tar ball](https://github.com/JeffreyRacine/R-Package-bkcde/tarball/main), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version:

```r
library(devtools); install_github('R-Package-bkcde', 'JeffreyRacine')
```

Note Windows users have to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X users have to first install [Xcode](https://apps.apple.com/us/app/xcode/id497799835) and the command line tools (in OS X 10.9 or higher, once you have Xcode installed, open a terminal and run xcode-select --install).
