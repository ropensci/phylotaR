#' @name genClRcrdBx
#' @title Generate Box of Cluster Records
#' @description Takes a list of ClRcrds, returns a box.
#' @param clstr_rcrds Cluster records
genClRcrdBx <- function(clstr_rcrds) {
  ids <- as.character(seq_along(clstr_rcrds)-1)
  names(clstr_rcrds) <- ids
  new('ClRcrdBx', ids=ids, cls=clstr_rcrds)
}

#' @name jnBxs
#' @title Join two cluster record boxes
#' @description Takes two ClRcrdBxs and joins them into
#' a single ClRcrdBx
#' @param bx_1 Cluster record box
#' @param bx_2 Cluster record box
jnBxs <- function(bx_1, bx_2) {
  genClRcrdBx(c(bx_1@cls, bx_2@cls))
}