#' @name fltr
#' @title Get all node IDs that will be processed
#' @description All nodes for which all children
#' contain < mxsqs sequences
fltr <- function(txid, phylt_nds, ps) {
  res <- vector()
  queue <- txid
  while(length(queue) > 0) {
    tmp_id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    nnonmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                              'n_gi_sub_nonmodel']
    nmodelsqs <- phylt_nds[match(tmp_id, phylt_nds[['ti']]),
                           'n_sp_model'] * ps[['mdlthrs']]
    nsqs <- nnonmodelsqs + nmodelsqs
    if(nsqs <= ps[['mxsqs']]) {
      res <- c(res, tmp_id)
    } else {
      info(lvl=2, ps=ps, "[", nsqs, " sqs] for [id ",
           tmp_id, "] ... searching descendants instead\n")
      queue <- c(queue, getDDFrmPhyltNds(tmp_id, phylt_nds))
    }
  }
  res
}

#' @name dwnld
#' @title Download sequences for txids
#' @description Look up and download all sequences for
#' given taxonomic IDs.
#' @param txids Taxonomic node IDs, numeric vector
#' @param phylt_nds PhyLoTa data.frame
#' @param txdct Taxonomic dictionary
#' @param verbose Verbose? T/F
#' @details Sequence downloads are cached.
dwnld <- function(txids, phylt_nds, txdct, ps) {
  # TODO: add overwrite arg
  sqcnt <- 0
  for(i in seq_along(txids)) {
    txid <- txids[i]
    info(lvl=1, ps=ps,
         "Downloading for [id ", txid, "]: [", i, "/",
         length(txids), "] ...")
    sqs <- getSqsByTxid(txid=txid, phylt_nds=phylt_nds,
                        ps=ps)
    if(length(sqs) > 0) {
      sqcnt <- sqcnt + length(sqs)
      sqs <- agmntSqRcrds(sqs=sqs, txdct=txdct)
      svSqs(wd=ps[['wd']], txid=txid, sqs=sqs)
    }
  }
  info(lvl=1, ps=ps, "Successfully downloaded [",
       sqcnt, " sqs] in total.")
}
