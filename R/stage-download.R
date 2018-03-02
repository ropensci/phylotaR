#' @name cldIdntfy
#' @title Get all node IDs that will be processed
#' @description All nodes for which all children
#' contain < mxsqs sequences
#' @param txdct TxDct
#' @param ps Parameters
cldIdntfy <- function(txdct, ps) {
  res <- vector()
  queue <- ps[['txid']]
  while(length(queue) > 0) {
    tmp_id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    sqcnt <- nSqs(txid=tmp_id, ps=ps)
    if(sqcnt <= ps[['mxsqs']]) {
      res <- c(res, tmp_id)
    } else {
      info(lvl=2, ps=ps, "[", nsqs, " sqs] for clade [id ",
           tmp_id, "] ... searching descendants instead\n")
      queue <- c(queue, getDDs(id=tmp_id, txdct=txdct))
    }
  }
  res
}

#' @name dwnld
#' @title Download sequences for txids
#' @description Look up and download all sequences for
#' given taxonomic IDs.
#' @param txids Taxonomic node IDs, numeric vector
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
    sqs <- getSqsByTxid(txid=txid, ps=ps)
    if(length(sqs) > 0) {
      sqcnt <- sqcnt + length(sqs)
      sqs <- agmntSqRcrds(sqs=sqs, txdct=txdct)
      svSqs(wd=ps[['wd']], txid=txid, sqs=sqs)
    }
  }
  info(lvl=1, ps=ps, "Successfully downloaded [",
       sqcnt, " sqs] in total.")
}
