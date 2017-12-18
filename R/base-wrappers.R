# wrappers for base functions
# required for testing as they need to be mocked
.system <- function(command, stdout, stderr) {
  base::system2(command, stdout, stderr)
}
.untar <- function(...) {
  utils::untar(...)
}