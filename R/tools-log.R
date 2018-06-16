# TODO: generate class for these?
#' @name info
#' @title Write info message to log
#' @description Inform a user via log.txt of pipeline progress.
#' @template lvl
#' @template ps
#' @param ... Message elements for concatenating
#' @family run-private
#' @return NULL
info <- function(lvl, ps, ...) {
  msg <- paste0(..., '\n')
  spcr <- paste0(rep('. ', lvl - 1), collapse = '')
  msg <- paste0(spcr, msg, collapse = '')
  .log(v = ps[['v']], wd = ps[['wd']], msg = msg)
}

#' @name error
#' @title Write error message to log
#' @description Inform a user if an error has occurred in log.txt,
#' halt pipeline.
#' @template ps
#' @param ... Message elements for concatenating
#' @family run-private
#' @return NULL
error <- function(ps, ...) {
  msg <- paste0(..., '\n')
  .log(v = FALSE, wd = ps[['wd']], msg = msg)
  stop(msg, call. = FALSE)
}

#' @name warn
#' @title Write warning message to log
#' @description Inform a user if a potential error has occurred in log.txt.
#' @template ps
#' @param ... Message elements for concatenating
#' @family run-private
#' @return NULL
warn <- function(ps, ...) {
  msg <- paste0('Warning: ', ..., '\n')
  .log(v = ps[['v']], wd = ps[['wd']], msg = msg)
  warning(paste0(msg, ' -- see log.txt'))
}

# hidden log function
.stgMsg <- function(ps, msg) {
  brdr <- paste0(rep('-', nchar(msg)), collapse = '')
  msg <- paste0(brdr, '\n', msg, '\n', brdr)
  info(ps = ps, lvl = 1, msg)
}
.log <- function(v, wd, msg) {
  if (v) {
    cat(msg)
  }
  if (!is.null(wd)) {
    lgfl <- file.path(wd, 'log.txt')
    cat(msg, file = lgfl, append = TRUE)
  }
}
