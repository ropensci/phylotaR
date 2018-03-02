#' @name getSqsByTxid
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a
#' taxonomic ID.
#' @param txid Taxonomic node ID, numeric
getSqsByTxid <- function(txid, txdct, ps, lvl=0) {
  # get subtree counts if that is smaller than ps[['mdlthrs']]
  subtree_count <- nSqs(txid, drct=FALSE, ps=ps)
  if(subtree_count <= ps[['mdlthrs']]) {
    info(lvl=2+lvl, ps=ps,
         '[', subtree_count, " sqs]. Downloading whole subtree.")
    sqs <- btchDwnld(txid=txid, drct=FALSE, ps=ps, lvl=lvl)
    return(sqs)
  }
  # 1st direct sqs from focal taxon, then from DDs
  sqs <- btchDwnld(txid=txid, drct=TRUE, ps=ps, lvl=lvl)
  for(dd in getDDs(id=txid, txdct=txdct)) {
    lvl <- lvl + 1
    info(lvl=2+lvl, ps=ps, "Downloading for child [id ", dd,"]")
    sqs <- c(sqs, getSqsByTxid(txid=dd, txdct=txdct,
                               ps=ps, lvl=lvl))
  }
  sqs
}

#' @name agmntSqRcrds
#' @title Augment sequence records list
#' @description Convert to sqsrcrds and add taxids
#' @param sqs List of SqRcrds
#' @param txdct Taxonomic Dictionary
agmntSqRcrds <- function(sqs, txdct) {
  # TODO: gen. txdct tools
  # txids are not downloaded as part of sequence, added here
  txdct_nms <- vapply(txdct@rcrds, function(x) x[['ScientificName']], '')
  txdct_ids <- vapply(txdct@rcrds, function(x) x[['TaxId']], '')
  sqs_nms <- vapply(sqs, function(x) x@orgnsm, '')
  sqs_ids <- txdct_ids[match(sqs_nms, txdct_nms)]
  for(i in seq_along(sqs)) {
    sqs[[i]]@txid <- as.character(sqs_ids[[i]])
  }
  genSqRcrdBx(sqs)
}
