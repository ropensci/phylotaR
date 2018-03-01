#' @name getSqsByTxid
#' @title Hierarchically get sequences for a txid
#' @description Looks up and downloads sequences for a
#' taxonomic ID.
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa data.frame
getSqsByTxid <- function(txid, phylt_nds, ps, lvl=0) {
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
  for(dd in getDDFrmPhyltNds(txid, phylt_nds)) {
    lvl <- lvl + 1
    info(lvl=2+lvl, ps=ps, "Downloading for child [id ", dd,"]")
    sqs <- c(sqs, getSqsByTxid(txid=dd, phylt_nds=phylt_nds,
                               ps=ps, lvl=lvl))
  }
  sqs
}

#' @name agmntRcrds
#' @title Augment sequence records list
#' @description Convert to sqsrcrds and add taxids
#' @param sqs List of SqRcrds
#' @param txdct Taxonomic Dictionary
agmntRcrds <- function(sqs, txdct) {
  # TODO: gen. txdct tools
  # txids are not downloaded as part of sequence, added here
  txdct_nms <- vapply(txdct, function(x) x[['ScientificName']], '')
  txdct_ids <- vapply(txdct, function(x) x[['TaxId']], '')
  sqs_nms <- vapply(sqs, function(x) x@orgnsm, '')
  sqs_ids <- txdct_ids[match(sqs_nms, txdct_nms)]
  for(i in seq_along(sqs)) {
    sqs[[i]]@txid <- as.character(sqs_ids[[i]])
  }
  genSqsRcrd(sqs)
}

#' @name getDDFrmPhyltNds
#' @title Get direct descendants from PhyLoTa nodes
#' @description Find next node IDs using the PhyLoTa
#' nodes data.frame.
#' @param txid parent
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details Returns vector of node IDs
getDDFrmPhyltNds <- function(txid, phylt_nds) {
  phylt_nds[phylt_nds[,'ti_anc']==txid,'ti']
}
