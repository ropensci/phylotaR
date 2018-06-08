#' @name cmdln
#' @title Run a command via terminal/command prompt
#' @description Provide the command and arguments as a vector.
#' Also can take a lgfl to which all stdout and stderr is written.
#' If lgfl is not provided, a list is returned of 'status', 'stdout'
#' and 'stderr'. Else only the status is returned - 1 success, 0
#' failed.
#' @param cmd Command to be run
#' @param args Vector of command arguments, each parameter and value
#' must be a separate element
#' @param lgfl File to which stdout/err will be written
#' @details Note, stdout/err are returned as 'raw'. Use rawToChar() to
#' convert to characters.
#' @family run-private
#' @return status, integer or character
cmdln <- function(cmd, args, lgfl=NULL) {
  if (!is.null(lgfl)) {
    # remove any filetype and replace with .log
    lgfl <- paste0(sub('\\.[^.]+$', '', lgfl), '.log')
    res <- try(sys::exec_wait(cmd = cmd, args = args,
                              std_out = lgfl, std_err = lgfl),
               silent = TRUE)
    if (inherits(res, 'try-error')) {
      cat(as.character(res), file = lgfl)
      res <- 1
    }
  } else {
    res <- try(sys::exec_internal(cmd = cmd, args = args),
               silent = TRUE)
    if (inherits(res, 'try-error')) {
      stderr <- charToRaw(as.character(res))
      res <- list('status' = 1, 'stderr' = stderr,
                  'stdout' = charToRaw(''))
    }
  }
  res
}

# wrappers for base functions
# required for testing as they need to be mocked
.untar <- function(...) {
  utils::untar(...)
}