#' @name getSqsByTxid
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a
#' taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa data.frame
#' @param verbose Verbose? T/F
#' @export
getSqsByTxid <- function(txid, phylt_nds, ps) {
  info(lvl=2, ps=ps, "Retrieving sequences for taxid [",
       txid, "]")
  # get subtree counts if that is smaller than ps[['mdlt']]
  subtree_count <- nSqs(txid, direct=FALSE, ps=ps)
  if(subtree_count <= ps[['mdlt']]) {
    info(lvl=3, ps=ps,
         '[', subtree_count, "] seqs for taxid [",
        txid, "], less than maximum of [", ps[['mdlt']],
        "] sequences. Retreiving sequences for whole subtree.")
    sqs <- dwnldFrmNCBI(txid=txid, direct=FALSE, ps=ps)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  sqs <- dwnldFrmNCBI(txid=txid, direct=TRUE, ps=ps)
  for(dd in getDDFrmPhyltNds(txid, phylt_nds)) {
    sqs <- c(sqs, getSqsByTxid(txid=dd, phylt_nds=phylt_nds,
                               ps=ps))
  }
  sqs
}

#' @name dwnld
#' @title Download sequences for txids
#' @description Look up and download all sequences for
#' given taxonomic IDs.
#' @param txids Taxonomic node IDs, numeric vector
#' @param phylt_nds PhyLoTa data.frame
#' @param ps[['mdlt']] Maximum number of sequences per txid
#' @param ps[['mxsql']] Maximum sequence length
#' @param verbose Verbose? T/F
#' @details Sequence downloads are cached.
#' @export
dwnld <- function(txids, phylt_nds, ps) {
  # TODO: add overwrite arg
  for(i in seq_along(txids)) {
    txid <- txids[i]
    info(lvl=2, ps=ps,
        "Downloading for taxid [", txid, "]: [", i, "/",
        length(txids), "]")
    sqs <- getSqsByTxid(txid=txid, phylt_nds=phylt_nds,
                        ps=ps)
    svSqs(wd=ps[['wd']], txid=txid, sqs=sqs)
    # Move this sqdf to separate function
    #sqdf <- do.call(rbind, lapply(sqs, as.data.frame))
  }
}

#' @name fltr
#' @title Get all node IDs that will be processed
#' @description All nodes for which all children
#' contain < mxsqs sequences
#' @export
fltr <- function(txid, phylt_nds, ps) {
  res <- vector()
  queue <- txid
  while(length(queue) > 0) {
    tmp_id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    # get number of sequences to determine if it
    # is manageable to calculate clusters
    info(lvl=3, ps=ps, "Counting species for taxon [",
        tmp_id, "]")
    nnonmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                               'n_gi_sub_nonmodel']
    nmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                            'n_sp_model'] * ps[['mdlt']]
    nsqs <- nnonmodelsqs + nmodelsqs
    info(lvl=3, ps=ps, "Number of sequences for taxon [",
        tmp_id, "]: [", nsqs, "]")
    # if sequence count is smaller than ps[['mxsqs']],
    # add it to the phylt_nds to process
    # otherwise look up direct descendents
    if(nsqs <= ps[['mxsqs']]) {
      info(lvl=3, ps=ps, "Will process taxon [",
           tmp_id, "]")
      res <- c(res, tmp_id)
    } else {
      info(lvl=3, ps=ps, "Too many seqs to blast for taxid [",
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
