
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
  prmtrs <- ldPrmtrs(wd)
  txid <- prmtrs[['txid']]
  mx_dscndnts <- prmtrs[['mx_dscndnts']]
  mx_sq_lngth <- prmtrs[['mx_sq_lngth']]
  mdl_thrshld <- prmtrs[['mdl_thrshld']]
  tdpth <- prmtrs[['tdpth']]
  tmout <- prmtrs[['tmout']]
  verbose <- prmtrs[['verbose']]
  # stage print
  msg <- paste0('Starting stage TAXISE: [', Sys.time(), ']')
  .stgMsg(v=verbose, wd=wd, msg=msg)
  # Run
  dwnldTD(wd=wd, tdpth=tdpth, verbose=verbose)
  tdobj <- genTDObj(wd=wd, verbose=verbose)
  info(lvl=1, v=verbose, wd=wd, 'Processing IDs ...')
  nid_sets <- getMngblIds(txid=txid,
                          td_nds=tdobj[['nds']],
                          mx_dscndnts=mx_dscndnts,
                          tmout=tmout,
                          verbose=verbose)
  info(lvl=1, v=verbose, wd=wd, 'Initiating PhyLoTa nodes ...')
  phylt_nds <- genPhylotaNds(wd=wd, nid_sets=nid_sets,
                             mx_sq_lngth=mx_sq_lngth,
                             mdl_thrshld=mdl_thrshld,
                             td_nds=tdobj[['nds']],
                             td_nms=tdobj[['nms']],
                             verbose=verbose)
  info(lvl=1, v=verbose, wd=wd, 'Writing out ...')
  svObj(wd=wd, obj=phylt_nds, nm='phylt_nds')
  writeTax(phylt_nds=phylt_nds, td_nms=tdobj[['nms']],
           fl=file.path(wd, paste0('dbfiles-taxonomy-',
                                   txid, '.tsv')),
           verbose=verbose)
  # stage print
  msg <- paste0('Completed stage TAXISE: [', Sys.time(), ']')
  .stgMsg(v=verbose, wd=wd, msg=msg)
}

#' @name runDownload
#' @title Download 
#' @description Download sequences
#' @export
runDownload <- function(wd) {
  # Get params
  prmtrs <- ldPrmtrs(wd)
  txid <- prmtrs[['txid']]
  mdl_thrshld <- prmtrs[['mdl_thrshld']]
  mx_blst_sqs <- prmtrs[['mx_blst_sqs']]
  mx_sq_lngth <- prmtrs[['mx_sq_lngth']]
  verbose <- prmtrs[['verbose']]
  # Get PhyLoTa nodes
  phylt_nds <- ldObj(wd=wd, nm='phylt_nds')
  # Filter
  fltrd_ids <- fltr(txid=txid, phylt_nds=phylt_nds,
                    mdl_thrshld=mdl_thrshld,
                    mx_blst_sqs=mx_blst_sqs,
                    verbose=verbose)
  # Download seqs
  dwnld(wd=wd, txids=fltrd_ids, phylt_nds=phylt_nds,
        mdl_thrshld=mdl_thrshld, mx_sq_lngth=mx_sq_lngth,
        verbose=verbose)
}

#' @name runClusters
#' @title Cluster
#' @description Identify sequence clusters
#' for an initial PhyLoTa cluster table
#' generated
#' @export
runClusters <- function(wd) {
  # Get params
  prmtrs <- ldPrmtrs(wd)
  txid <- prmtrs[['txid']]
  mdl_thrshld <- prmtrs[['mdl_thrshld']]
  mx_blst_sqs <- prmtrs[['mx_blst_sqs']]
  mx_sq_lngth <- prmtrs[['mx_sq_lngth']]
  verbose <- prmtrs[['verbose']]
  # Get PhyLoTa nodes
  phylt_nds <- ldObj(wd=wd, nm='phylt_nds')
  # generate clusters
  calcClstrs(wd=wd, txid=txid, phylt_nds=phylt_nds,
             verbose=verbose)
}

#' @name runAlign
#' @title Align clusters
#' @description Run external alignment software
#' on identified clusters.
#' @param wd Working directory
#' @details Object will be cached.
#' @export
runAlign <- function(wd) {
  # TODO
}