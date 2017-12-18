# wrappers for base functions
# required for testing as they need to be mocked
.system <- function(command, args, stdout, stderr) {
  base::system2(command, args, stdout, stderr)
}
.untar <- function(...) {
  utils::untar(...)
}