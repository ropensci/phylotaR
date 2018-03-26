# GEN
#' @name genTxTr
#' @title Generate taxonomic tree
#' @description Generate a taxonomic tree for
#' easy look up of taxonomic parents and descendants.
#' @return TreeMan
#' @param prinds Vector of integers indicating preceding node.
#' @param trids Vector of taxonomic IDs
#' @param root ID of root taxon
#' @importClassesFrom treeman TreeMan
genTxTr <- function(prinds, trids, root) {
  .add <- function(i) {
    nd <- vector("list", length=4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- trids[i]
    nd[['prid']] <- trids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd[['spn']] <- 1
    nd
  }
  nonroot_i <- trids != root
  nnds <- length(prinds)
  tinds <- which(!1:nnds %in% prinds)
  ptnds_pool <- prinds[nonroot_i]
  ptids <- trids[nonroot_i]
  ndlst <- lapply(1:nnds, .add)
  names(ndlst) <- trids
  tree <- new('TreeMan', ndlst=ndlst, root=root,
              wtxnyms=FALSE, ndmtrx=NULL,
              prinds=prinds, tinds=tinds)
  treeman::checkNdlst(ndlst, root)
  treeman::fastCheckTreeMan
  treeman::updateSlts(tree)
}

# CONVERSION
#' @name trid2Txid
#' @title Get taxonomic ID from tree ID
#' @description Return taxonomic ID for taxonomic
#' tree ID in TxDct.
#' @return Character or vector
#' @param id txid
#' @param txdct TxDct
trid2Txid <- function(id, txdct) {
  txdct@txids[match(id, txdct@indx)]
}

#' @name txid2Trid
#' @title Get tree ID from taxonomic ID
#' @description Return taxonomic tree ID from taxonomic
#' ID in TxDct.
#' @return Character or vector
#' @param id trid
#' @param txdct TxDct
txid2Trid <- function(id, txdct) {
  as.character(txdct@indx[match(id, txdct@txids)])
}

# GET
#' @name getSngltns
#' @title Get equivalent singleton IDs
#' @description For a given taxonomic ID, return all the
#' equivalent singleton IDs.
#' @return Character vector
#' @param txid txid
#' @param txdct TxDct
#' @details phylotaR drops all singleton nodes from
#' the NCBI taxonomy upon initiation. To retrieve all
#' the singleton IDs for an ID use this function.
#' Singletons are defined as nodes in the hierarchy
#' that do not split. 
getSngltns <- function(txid, txdct) {
  trid <- as.character(txdct@indx[match(txid, txdct@txids)])
  sngltns <- txdct@txids[trid == txdct@indx]
  sngltns[sngltns != txid]
}

#' @name getRnk
#' @title Get rank
#' @description Look-up taxonomic rank from TxDct.
#' @return Character
#' @param id txid
#' @param txdct TxDct
getRnk <- function(id, txdct) {
  txdct@rcrds[[id]]@rnk
}

#' @name getADs
#' @title Get all descendants
#' @description Look-up all taxonomic descendants of a node from TxDct.
#' @return Vector
#' @param id txid
#' @param txdct TxDct
getADs <- function(id, txdct) {
  trid <- txid2Trid(id=id, txdct=txdct)
  tr_ptids <- treeman::getNdPtids(tree=txdct@txtr, id=trid)
  trid2Txid(id=tr_ptids, txdct=txdct)
}

#' @name getDDs
#' @title Get direct descendants
#' @description Look-up direct taxonomic descendants of a node from TxDct.
#' @return Vector
#' @param id txid
#' @param txdct TxDct
getDDs <- function(id, txdct) {
  trid <- txid2Trid(id=id, txdct=txdct)
  tr_ptids <- treeman::getNdSlt(tree=txdct@txtr,
                                slt_nm='ptid', id=trid)
  trid2Txid(id=tr_ptids, txdct=txdct)
}

#' @name getPrnt
#' @title Get taxonomic parent
#' @description Look-up MRCA of id(s) from TxDct
#' @return Character
#' @param id txid(s)
#' @param txdct TxDct
getPrnt <- function(id, txdct) {
  trid <- txid2Trid(id=id, txdct=txdct)
  if(length(trid) > 1) {
    res <- treeman::getPrnt(tree=txdct@txtr,
                            ids=unique(trid))
  } else {
    res <- treeman::getNdSlt(tree=txdct@txtr,
                             slt_nm='prid', id=trid)
  }
  trid2Txid(id=res, txdct=txdct)
}
