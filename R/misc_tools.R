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

# lvl - number of indentations, indicating code depth
# v - verbose, T/F
# wd - working directory for accessing log file
# ... message parts
.log <- function(lvl, v, wd, ...) {
  msg <- paste0(..., '\n')
  spcr <- paste0(rep('... ', lvl-1), collapse='')
  msg <- paste0(spcr, msg, collapse='')
  if(v) {
    cat(msg)
  }
  if(!is.null(wd)) {
    lgfl <- file.path(wd, 'log.txt')
    cat(msg, file=lgfl, append=TRUE)
  }
}

# start and end stage messages
.stgLog <- function(v, wd, msg) {
  brdr <- paste0(rep('-', nchar(msg)), collapse='')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  .log(v=v, wd=wd, lvl=1, msg)
}