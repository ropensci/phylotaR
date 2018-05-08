#' @name archive_cluster_generate
#' @title Generate ArchiveCluster container class
#' @description Takes a list of RecordCluster classes, returns an ArchiveCluster.
#' @param clstr_rcrds list of RecordCluster classes
#' @return ArchiveCluster
#' @noRd
archive_cluster_generate <- function(clstr_rcrds) {
  ids <- as.character(seq_along(clstr_rcrds) - 1)
  names(clstr_rcrds) <- ids
  new('ArchiveCluster', ids = ids, cls = clstr_rcrds)
}

#' @name archive_cluster_join
#' @title Join two ArchiveCluster classes
#' @description Take two ArchiveCluster classes and join them into
#' a single ArchiveCluster
#' @param ac_1 ArchiveCluster
#' @param ac_2 ArchiveCluster
#' @return ArchiveCluster
#' @noRd
archive_cluster_join <- function(ac_1, ac_2) {
  archive_cluster_generate(c(ac_1@cls, ac_2@cls))
}