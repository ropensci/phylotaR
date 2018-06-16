#' @name download_run
#' @title Run download stage
#' @description Run the second stage of phylotaR, download. This stage
#' downloads sequences for all nodes with sequence numbers less than
#' mxsqs. It hierarchically traverses the taxonomy for each node and
#' downloads direct and subtree sequences for all descendants.
#' @param wd Working directory
#' @export
#' @example examples/download_run.R
#' @return NULL
download_run <- function(wd) {
  ps <- parameters_load(wd)
  msg <- paste0('Starting stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  info(lvl = 1, ps = ps, 'Identifying suitable clades ...')
  txdct <- obj_load(wd = wd, nm = 'txdct')
  clds_ids <- clade_select(txdct = txdct, ps = ps)
  info(lvl = 1, ps = ps, 'Identified [', length(clds_ids),
       '] suitable clades.')
  info(lvl = 1, ps = ps, 'Downloading hierarchically ...')
  seq_download(txids = clds_ids, txdct = txdct, ps = ps)
  msg <- paste0('Completed stage DOWNLOAD: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name clade_select
#' @title Get all node IDs that will be processed
#' @description All nodes with less than maximum number
#' of nodes and sequences.
#' @param txdct TxDct
#' @template ps
#' @return vector of txids
#' @family run-private
clade_select <- function(txdct, ps) {
  res <- vector()
  queue <- ps[['txid']]
  while (length(queue) > 0) {
    tmp_id <- head(queue, 1)
    queue <- tail(queue, length(queue) - 1)
    sqcnt <- sqs_count(txid = tmp_id, ps = ps)
    ndcnt <- txnds_count(txid = tmp_id, ps = ps)
    mx_pssbl_sqcnt <- ps[['mdlthrs']] * ndcnt
    sqcnt <- ifelse(sqcnt > mx_pssbl_sqcnt, mx_pssbl_sqcnt, sqcnt)
    if (sqcnt <= ps[['mxsqs']] & ndcnt <= ps[['mxnds']]) {
      res <- c(res, tmp_id)
    } else {
      info(lvl = 2, ps = ps, "[", sqcnt, " sqs] and [", ndcnt,
           " nds] for clade [id ", tmp_id,
           "] - searching descendants instead ...")
      queue <- c(queue, descendants_get(id = as.character(tmp_id),
                                        direct = TRUE, txdct = txdct))
    }
  }
  res
}

#' @name seq_download
#' @title Download sequences for txids
#' @description Look up and download all sequences for given
#' taxonomic IDs.
#' @param txids Taxonomic node IDs, numeric vector
#' @param txdct Taxonomic dictionary
#' @template ps
#' @details Sequence downloads are cached.
#' @return NULL
#' @family run-private
seq_download <- function(txids, txdct, ps) {
  # TODO: add overwrite arg
  sqcnt <- 0
  for (i in seq_along(txids)) {
    txid <- txids[i]
    sqfl <- file.path(ps[['wd']], 'cache', 'sqs',
                      paste0(txid, '.RData'))
    if (file.exists(sqfl)) {
      next
    }
    info(lvl = 1, ps = ps, "Working on parent [id ", txid, "]: [", i,
         "/", length(txids), "] ...")
    sqs <- hierarchic_download(txid = txid, txdct = txdct, ps = ps)
    if (length(sqs) > 0) {
      sqcnt <- sqcnt + length(sqs)
      sqs <- seqrec_augment(sqs = sqs, txdct = txdct)
      sqs_save(wd = ps[['wd']], txid = txid, sqs = sqs)
    }
  }
  info(lvl = 1, ps = ps, "Successfully downloaded [", sqcnt,
       " sqs] in total.")
}
