.onAttach <- function (lib, pkg) {
	packageStartupMessage("package bkcde version 1.47")
}

## Declare globals used in data.frame operations to silence R CMD check NOTES
if(getRversion() >= "2.15.1") utils::globalVariables(
  c("error", "optim_elapsed", "convergence", "elapsed_fit", "elapsed_total")
)
