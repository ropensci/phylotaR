
#' @name runTaxise
#' @title Run taxise stage
#' @description TODO
#' @param wd Working directory
#' @details Object will be cached.
#' @export
runTaxise <- function(wd) {
  # TODO: allow a user to have their own taxids and/or tax tree
  ps <- ldPrmtrs(wd)
  msg <- paste0('Starting stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  info(lvl=1, ps=ps, 'Searching taxonomic IDs ...')
  txids <- getTxids(ps=ps)
  info(lvl=1, ps=ps, 'Downloading taxonomic records ...')
  rcrds <- dwnldTxRcrds(txids=txids, ps=ps)
  info(lvl=1, ps=ps, 'Generating taxonomic dictionary ...')
  txdct <- genTxDct(rcrds=rcrds, txids=txids)
  svObj(wd=wd, obj=txdct, nm='txdct')
  msg <- paste0('Completed stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runDownload
#' @title Download 
#' @description Download sequences
#' @param wd Working directory
#' @export
runDownload <- function(wd) {
  ps <- ldPrmtrs(wd)
  msg <- paste0('Starting stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  info(lvl=1, ps=ps, 'Identifying suitable clades ...')
  txdct <- ldObj(wd=wd, nm='txdct')
  clds_ids <- cldIdntfy(txdct=txdct, ps=ps)
  info(lvl=1, ps=ps, 'Identified [', length(clds_ids),
       '] suitable clades.')
  info(lvl=1, ps=ps, 'Downloading hierarchically ...')
  dwnldSqRcrds(txids=clds_ids, txdct=txdct, ps=ps)
  msg <- paste0('Completed stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runClusters
#' @title Cluster
#' @description Identify sequence clusters
#' for an initial PhyLoTa cluster table
#' generated
#' @param wd Working directory
#' @export
runClusters <- function(wd) {
  ps <- ldPrmtrs(wd)
  msg <- paste0('Starting stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  txdct <- ldObj(wd=wd, nm='txdct')
  calcClstrs(ps=ps, txdct=txdct)
  msg <- paste0('Completed stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runClusters2
#' @title Cluster
#' @description Cluster the clusters
#' @param wd Working directory
#' @export
runClusters2 <- function(wd) {
  ps <- ldPrmtrs(wd)
  # stage print
  msg <- paste0('Starting stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  # generate clusters
  clstrClstrs(ps=ps)
  # stage print
  msg <- paste0('Completed stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}
