
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
  tmout <- prmtrs[['tmout']]
  verbose <- prmtrs[['verbose']]
  # Run
  tdobj <- genTDObj(wd)
  cat('Processing IDs....\n')
  nid_sets <- getMngblIds(txid=txid,
                          td_nds=tdobj[['nds']],
                          mx_dscndnts=mx_dscndnts,
                          tmout=tmout,
                          verbose=verbose)
  cat('Initiating PhyLoTa nodes....\n')
  phylt_nds <- genPhylotaNds(nid_sets=nid_sets,
                             mx_sq_lngth=mx_sq_lngth,
                             mdl_thrshld=mdl_thrshld,
                             td_nds=tdobj[['nds']],
                             td_nms=tdobj[['nms']],
                             verbose=verbose)
  cat('Writing out....\n')
  svObj(wd=wd, obj=phylt_nds, nm='phylt_nds')
  writeTax(phylt_nds=phylt_nds, td_nms=tdobj[['nms']],
           fl=file.path(wd, paste0('dbfiles-taxonomy-',
                                   txid, '.tsv')),
           verbose=verbose)
  cat('Done.\n')
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
