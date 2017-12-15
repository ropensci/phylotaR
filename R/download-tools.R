#' @name getSqsByTxid
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a
#' taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa data.frame
#' @param mx_len Maximum sequence length
#' @param mdl_threshold Maximum sequences per species
#' @param verbose Verbose? T/F
#' @export
# TODO: make uniform mx_len
getSqsByTxid <- function(wd, txid, phylt_nds, mx_len=25000,
                         mdl_thrshld=10000, verbose=FALSE) {
  info(lvl=2, v=verbose, wd=wd, "Retrieving sequences for taxid [",
       txid, "]")
  # get subtree counts if that is smaller than mdl_thrshld
  subtree_count <- nSqs(txid, direct=FALSE,
                        mx_len=mx_len, verbose=verbose)
  if(subtree_count <= mdl_thrshld) {
    info(lvl=3, v=verbose, wd=wd,
         '[', subtree_count, "] seqs for taxid [",
        txid, "], less than maximum of [", mdl_thrshld,
        "] sequences. Retreiving sequences for whole subtree.")
    sqs <- dwnldFrmNCBI(txid=txid, direct=FALSE, mx_lngth=mx_len,
                        mx_sqs=mdl_thrshld, verbose=verbose)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  sqs <- dwnldFrmNCBI(txid=txid, direct=TRUE, mx_lngth=mx_len,
                      mx_sqs=mdl_thrshld, verbose=verbose)
  for(dd in getDDFrmPhyltNds(txid, phylt_nds)) {
    sqs <- c(sqs, getSqsByTxid(txid=dd, phylt_nds=phylt_nds,
                               mx_len=mx_len,
                               mdl_thrshld=mdl_thrshld,
                               verbose=verbose))
  }
  sqs
}

#' @name dwnld
#' @title Download sequences for txids
#' @description Look up and download all sequences for
#' given taxonomic IDs.
#' @param txids Taxonomic node IDs, numeric vector
#' @param phylt_nds PhyLoTa data.frame
#' @param mdl_thrshld Maximum number of sequences per txid
#' @param mx_sq_lngth Maximum sequence length
#' @param verbose Verbose? T/F
#' @details Sequence downloads are cached.
#' @export
dwnld <- function(wd, txids, phylt_nds, mdl_thrshld,
                  mx_sq_lngth, verbose) {
  # TODO: add overwrite arg
  for(i in seq_along(txids)) {
    txid <- txids[i]
    info(lvl=2, v=verbose, wd=wd,
        "Downloading for taxid [", txid, "]: [", i, "/",
        length(txids), "]")
    sqs <- getSqsByTxid(txid=txid, phylt_nds=phylt_nds,
                        mx_len=mx_sq_lngth,
                        mdl_thrshld=mdl_thrshld,
                        verbose=verbose)
    svSqs(wd=wd, txid=txid, sqs=sqs)
    # Move this sqdf to separate function
    #sqdf <- do.call(rbind, lapply(sqs, as.data.frame))
  }
}

#' @name fltr
#' @title Get all node IDs that will be processed
#' @description All nodes for which all children
#' contain < mx_blst_sqs sequences
#' @export
fltr <- function(wd, txid, phylt_nds, mdl_thrshld,
                 mx_blst_sqs, verbose=FALSE) {
  res <- vector()
  queue <- txid
  while(length(queue) > 0) {
    tmp_id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    # get number of sequences to determine if it
    # is manageable to calculate clusters
    info(lvl=3, v=verbose, wd=wd, "Counting species for taxon [",
        tmp_id, "]")
    nnonmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                               'n_gi_sub_nonmodel']
    nmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                            'n_sp_model'] * mdl_thrshld
    nsqs <- nnonmodelsqs + nmodelsqs
    info(lvl=3, v=verbose, wd=wd, "Number of sequences for taxon [",
        tmp_id, "]: [", nsqs, "]")
    # if sequence count is smaller than mx_blst_sqs,
    # add it to the phylt_nds to process
    # otherwise look up direct descendents
    if(nsqs <= mx_blst_sqs) {
      info(lvl=3, v=verbose, wd=wd, "Will process taxon [",
           tmp_id, "]")
      res <- c(res, tmp_id)
    } else {
      info(lvl=3, v=verbose, wd=wd, "Too many seqs to blast for taxid [",
          tmp_id, "] ... looking up direct descendents\n")
      queue <- c(queue, getDDFrmPhyltNds(tmp_id, phylt_nds))
    }
  }
  res
}

#' @name getDDFrmPhyltNds
#' @title Get direct descendants from PhyLoTa nodes
#' @description Find next node IDs using the PhyLoTa
#' nodes data.frame.
#' @param txid parent
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details Returns vector of node IDs
#' @export
getDDFrmPhyltNds <- function(txid, phylt_nds) {
  phylt_nds[phylt_nds[,'ti_anc']==txid,'ti']
}
