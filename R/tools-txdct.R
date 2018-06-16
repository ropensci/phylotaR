#' @name taxtree_gen
#' @title Generate taxonomic tree
#' @description Generate a taxonomic tree for
#' easy look up of taxonomic parents and descendants.
#' @return TreeMan
#' @param prinds Vector of integers indicating preceding node.
#' @param ids Vector of taxonomic IDs
#' @param root ID of root taxon
#' @template ps
#' @importClassesFrom treeman TreeMan
#' @family run-private
#' @return TreeMan class
taxtree_gen <- function(prinds, ids, root, ps) {
  .add <- function(i) {
    nd <- vector("list", length = 4)
    names(nd) <- c('id', 'ptid', 'prid', 'spn')
    nd[['id']] <- ids[i]
    nd[['prid']] <- ids[prinds[i]]
    nd[['ptid']] <- ptids[ptnds_pool == i]
    nd[['spn']] <- 1
    nd
  }
  nonroot_i <- ids != root
  nnds <- length(prinds)
  tinds <- which(!1:nnds %in% prinds)
  ptnds_pool <- prinds[nonroot_i]
  ptids <- ids[nonroot_i]
  ndlst <- lapply(1:nnds, .add)
  names(ndlst) <- ids
  tree <- treeman::twoer()
  tree@ndlst <- ndlst
  tree@root <- root
  tree@wtxnyms <- FALSE
  tree@ndmtrx <- NULL
  tree@prinds <- prinds
  tree@tinds <- tinds
  tree <- treeman::updateSlts(tree)
  if(!treeman::checkNdlst(tree@ndlst, tree@root)) {
    error(ps=ps, 'Invalid taxonomy')
  }
  tree
}

#' @name rank_get
#' @title Get rank
#' @description Look-up taxonomic rank from dictionary.
#' @return character
#' @param txid txid
#' @param txdct TaxDict
#' @family run-private
rank_get <- function(txid, txdct) {
  txdct@recs[[txid]]@rnk
}

#' @name descendants_get
#' @title Get descendants
#' @description Look-up either direct or all taxonomic descendants of
#' a node from taxonomic dictionary.
#' @return vector
#' @param id txid
#' @param txdct TaxDict
#' @param direct T/F, return only direct descendants?
#' @family run-private
descendants_get <- function(id, txdct, direct=FALSE) {
  if (direct) {
    ptids <- treeman::getNdSlt(tree = txdct@txtr, slt_nm = 'ptid',
                               id = id)
  } else {
    ptids <- treeman::getNdPtids(tree = txdct@txtr, id = id)
  }
  ptids
}

#' @name parent_get
#' @title Get taxonomic parent
#' @description Look-up MRCA of taxonomic id(s) from taxonomic
#' dictionary
#' @return Character
#' @param id txid(s)
#' @param txdct TaxDict
#' @family run-private
parent_get <- function(id, txdct) {
  if (length(id) > 1) {
    res <- treeman::getPrnt(tree = txdct@txtr,
                            ids = unique(id))
  } else {
    res <- treeman::getNdSlt(tree = txdct@txtr,
                             slt_nm = 'prid', id = id)
  }
  res
}
