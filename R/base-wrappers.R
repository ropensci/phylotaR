# wrappers for base functions
# required for testing as they need to be mocked
.system <- function(command, args, fl=NULL) {
  if(!is.null(fl)) {
    stdfl <- paste0(sub('\\.fa', '', fl), '.log')
  } else {
    stdfl <- FALSE
  }
  sys::exec_wait(cmd=command, args=args, std_out=stdfl,
                 std_err=stdfl)
}
.untar <- function(...) {
  utils::untar(...)
}