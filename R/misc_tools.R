# Hidden functions

.evlTmLmt <- function(expr, cpu=Inf, elapsed=Inf) {
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent=TRUE)
  setTimeLimit()
  if(inherits(y, "try-error")) NULL else y
}

.cp <- function(v, ...) {
  # TODO: add levels arg
  # Custom print
  if(v) {
    cat(..., '\n', sep='')
  }
}