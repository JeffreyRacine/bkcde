# bkcde

This repository contains the R package `bkcde' written and maintained by Jeffrey S. Racine (racinej@mcmaster.ca)

## Installation

You can install this package by downloading the [zip ball](https://github.com/JeffreyRacine/bkcde/zipball/main) or [tar ball](https://github.com/JeffreyRacine/bkcde/tarball/main), decompress and run `R CMD INSTALL` on it, or install then use the **devtools** package to install the development version, as follows:

```r
install.packages(c("MCPAN","pbmcapply","robustbase"),repos="https://cran.wu.ac.at/")
library(devtools)
install_github('JeffreyRacine/bkcde')
```

Note Windows users may have to first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), while OS X users may have to first install [Xcode](https://apps.apple.com/us/app/xcode/id497799835) and the command line tools (in OS X 10.9 or higher, once you have Xcode installed, open a terminal and run xcode-select --install). If you are a regular R user you can probably skip this step.


