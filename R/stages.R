
#' @name runTaxise
#' @title Extract and filter NCBI taxonomy
#' @description Look up NCBI taxonomy and generate
#' a list of descendent taxonomic nodes based on
#' provided root \code{txid}. Descendent nodes
#' are written to disk in PhyLoTa table format.
#' This is a necessary step before clustering can
#' be performed.
#' @param wd Working directory
#' @details Object will be cached.
#' @export
runTaxise <- function(wd) {
  # Get params
  ps <- ldPrmtrs(wd)
  # stage print
  msg <- paste0('Starting stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  # Run
  dwnldTD(ps=ps)
  tdobj <- genTDObj(ps=ps)
  info(lvl=1, ps=ps, 'Processing IDs ...')
  nid_sets <- getMngblIds(txid=ps[['txid']], td_nds=tdobj[['nds']],
                          ps=ps)
  info(lvl=1, ps=ps, 'Initiating PhyLoTa nodes ...')
  phylt_nds <- genPhylotaNds(nid_sets=nid_sets,
                             td_nds=tdobj[['nds']],
                             td_nms=tdobj[['nms']],
                             ps=ps)
  info(lvl=1, ps=ps, 'Writing out ...')
  svObj(wd=wd, obj=phylt_nds, nm='phylt_nds')
  # TODO: make this optional?
  info(lvl=1, ps=ps, 'Generating taxonomic dictionary ...')
  txdct <- genTxdct(phylt_nds=phylt_nds, ps=ps)
  svObj(wd=wd, obj=txdct, nm='txdct')
  # stage print
  msg <- paste0('Completed stage TAXISE: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runDownload
#' @title Download 
#' @description Download sequences
#' @export
runDownload <- function(wd) {
  # Get params
  ps <- ldPrmtrs(wd)
  # stage print
  msg <- paste0('Starting stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  # Get PhyLoTa nodes
  phylt_nds <- ldObj(wd=wd, nm='phylt_nds')
  info(lvl=1, ps=ps, 'Filtering ...')
  fltrd_ids <- fltr(txid=ps[['txid']], phylt_nds=phylt_nds, ps=ps)
  info(lvl=1, ps=ps, 'Downloading ...')
  dwnld(txids=ps[['txid']], phylt_nds=phylt_nds, ps=ps)
  # stage print
  msg <- paste0('Completed stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runClusters
#' @title Cluster
#' @description Identify sequence clusters
#' for an initial PhyLoTa cluster table
#' generated
#' @export
runClusters <- function(wd) {
  ps <- ldPrmtrs(wd)
  # stage print
  msg <- paste0('Starting stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
  # Get PhyLoTa nodes
  phylt_nds <- ldObj(wd=wd, nm='phylt_nds')
  # generate clusters
  calcClstrs(txid=ps[['txid']], phylt_nds=phylt_nds, ps=ps)
  # stage print
  msg <- paste0('Completed stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps=ps, msg=msg)
}

#' @name runClusters2
#' @title Cluster
#' @description Cluster the clusters
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