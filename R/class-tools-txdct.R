# GEN
#' @name genTxTr
#' @title Generate taxonomic tree
#' @description Generate a taxonomic tree for
#' easy look up of taxonomic parents and descendants.
#' @return TreeMan
#' @param prinds Vector of integers indicating preceding node.
#' @param txids Vector of taxonomic IDs
#' @param root ID of root taxon
#' @importClassesFrom treeman TreeMan
genTxTr <- function(prinds, txids, root) {
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- txids[i]
    nd[['prid']] <- txids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd[['spn']] <- 1
    nd
  }
  nonroot_i <- txids != root
  nnds <- length(prinds)
  tinds <- which(!1:nnds %in% prinds)
  ptnds_pool <- prinds[nonroot_i]
  ptids <- txids[nonroot_i]
  ndlst <- lapply(1:nnds, .add)
  names(ndlst) <- txids
  tree <- new('TreeMan', ndlst=ndlst, root=root,
              wtxnyms=FALSE, ndmtrx=NULL,
              prinds=prinds, tinds=tinds)
  treeman::updateSlts(tree)
}

# GET
#' @name getRnk
#' @title Get rank
#' @description Look-up taxonomic rank from TxDct.
#' @return Character
#' @param id txid
#' @param txdct TxDct
getRnk <- function(id, txdct) {
  txdct@rcrds[[id]][['Rank']]
}

#' @name getADs
#' @title Get all descendants
#' @description Look-up all taxonomic descendants of a node from TxDct.
#' @return Vector
#' @param id txid
#' @param txdct TxDct
getADs <- function(id, txdct) {
  treeman::getNdPtids(tree=txdct@txtr, id=id)
}

#' @name getDDs
#' @title Get direct descendants
#' @description Look-up direct taxonomic descendants of a node from TxDct.
#' @return Vector
#' @param id txid
#' @param txdct TxDct
getDDs <- function(id, txdct) {
  treeman::getNdSlt(tree=txdct@txtr, slt_nm='ptid', id=id)
}

#' @name getPrnt
#' @title Get taxonomic parent
#' @description Look-up MRCA of id(s) from TxDct
#' @return Character
#' @param id txid(s)
#' @param txdct TxDct
getPrnt <- function(id, txdct) {
  if(length(id) > 1) {
    return(treeman::getPrnt(tree=txdct@txtr, ids=id))
  }
  treeman::getNdSlt(tree=txdct@txtr, slt_nm='prid', id=id)
}


#' @name getIDFrmTxdct
#' @title Get taxonomic IDs by rank
#' @description
#' @param txdct Taxonomic dictionary
#' @param id Taxon IDs
#' @param ret Return ID or Name?
#' @param rank Rank of output
getIDFrmTxdct <- function(txdct, id, ret=c('TaxId', 'ScientificName'),
                          rank=c("superkingdom", "kingdom", "phylum",
                                 "subphylum", "class", "superorder",
                                 "order", "suborder", "infraorder",
                                 "parvorder", "family", "genus",
                                 "species", "subspecies")) {
  calc <- function(id) {
    rcrd <- txdct[[id]]
    rnks <- sapply(rcrd[['LineageEx']],
                   function(x) x[['Rank']])
    ids <- sapply(rcrd[['LineageEx']],
                  function(x) x[[ret]])
    mtch <- which(rank == rnks)
    if(length(mtch) == 0) {
      # if no match, use lowest available rank
      mtch <- length(rnks)
    }
    ids[[mtch[[1]]]]
  }
  rank <- match.arg(rank)
  ret <- match.arg(ret)
  sapply(as.character(id), calc)
}