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
#' @param txid Root taxonomic NCBI ID
#' @details Object will be cached.
#' @export
# @Hannes: this is functional eqv to create.nodes
genTxNds <- function(wd, txid) {
  # Get params
  prmtrs <- ldPrmtrs(wd)
  mx_dscndnts <- prmtrs[['mx_dscndnts']]
  tmout <- prmtrs[['tmout']]
  verbose <- prmtrs[['verbose']]
  # Run
  tdobj <- genTDObj(wd)
  nd_ids <- getMngblIds(txid=txid,
                        td_nds=tdobj[['nds']],
                        mx_dscndnts=mx_dscndnts,
                        tmout=tmout,
                        verbose=verbose)
}

