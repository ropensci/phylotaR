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
#' @param ps Paramters
#' @details Note, stdout/err are returned as 'raw'. Use rawToChar() to
#' convert to characters.
#' @family run-private
#' @return status, integer or character
cmdln <- function(cmd, args, ps, lgfl = NULL) {
  if (!is.null(lgfl)) {
    # remove any filetype and replace with .log
    lgfl <- paste0(sub("\\.[^.]+$", "", lgfl), ".log")
    res <- sys_exec(cmd = cmd, args = args, ps = ps, lgfl = lgfl)
    if (inherits(res, "try-error")) {
      cat(as.character(res), file = lgfl)
      res <- 1
    }
  } else {
    res <- sys_exec(cmd = cmd, args = args, ps = ps)
    if (inherits(res, "try-error")) {
      stderr <- charToRaw(as.character(res))
      res <- list(
        "status" = 1, "stderr" = stderr,
        "stdout" = charToRaw("")
      )
    }
  }
  res
}

# wrappers for base functions
# required for testing as they need to be mocked
.untar <- function(...) {
  utils::untar(...)
}

# 9 Sept 2019 outsider integration
repo <- "dombennett/om..blast"
sys_exec <- function(cmd = cmd, args = args, ps, lgfl = NULL) {
  if (ps[["outsider"]]) {
    if (!outsider_call("is_module_installed")(repo = repo)) {
      message("Installing BLAST with outsider")
      outsider_call("module_install")(
        repo = repo, service = "github", tag = "latest",
        force = TRUE
      )
    }
    cmd <- outsider_call("module_import")(fname = cmd, repo = repo)
    if (!is.null(lgfl)) {
      outsider_call("verbosity_set")(show_program = lgfl)
    }
    res <- try(cmd(arglist = args), silent = TRUE)
    if (is.logical(res) && res) {
      res <- 0
    }
    if (is.null(lgfl)) {
      res <- list("status" = res, "stderr" = raw(), "stdout" = raw())
    }
  } else {
    if (!is.null(lgfl)) {
      res <- try(sys::exec_wait(
        cmd = cmd, args = args, std_out = lgfl,
        std_err = lgfl
      ), silent = TRUE)
    } else {
      res <- try(sys::exec_internal(cmd = cmd, args = args), silent = TRUE)
    }
  }
  res
}

outsider_call <- function(.f_string) {
  stopifnot(is.character(.f_string))
  cfun <- paste0("outsider::", .f_string)
  eval(parse(text = cfun))
}
