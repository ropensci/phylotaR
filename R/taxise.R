# Demo-ing the restructuring of the pipeline
# for better modularity, testing and package dev
# This file will hold the stage functions required
# for interacting with the NCBI taxonomy

#' @name genTxNds
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
# @Hannes: this is functional eqv to create.nodes
genTxNds <- function(wd) {
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
  writeTax(phylt_nds=phylt_nds, td_nms=tdobj[['nms']],
           fl=file.path(wd, paste0('dbfiles-taxonomy-', txid, '.tsv')),
           verbose=verbose)
  cat('Done.\n')
}

