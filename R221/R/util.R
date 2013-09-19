
## Get fractions of a second of the current time
fracSec <- function() {
  as.numeric(format(Sys.time(), "%OS")) * 10^5
}

## Run all tests of the package
runTests <- function() {
  require(testthat)
  tdir <- system.file("tests", package="R221")
  test_dir(tdir)
}
