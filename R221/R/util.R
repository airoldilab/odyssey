
## Get fractions of a second of the current time
fracSec <- function() {
  as.numeric(format(Sys.time(), "%OS")) * 10^5
}

