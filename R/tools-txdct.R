#' @name taxtree_gen
#' @title Generate taxonomic tree
#' @description Generate a taxonomic tree for
#' easy look up of taxonomic parents and descendants.
#' @return TreeMan
#' @param prinds Vector of integers indicating preceding node.
#' @param trids Vector of taxonomic IDs
#' @param root ID of root taxon
#' @importClassesFrom treeman TreeMan
#' @noRd
#' @return TreeMan class
taxtree_gen <- function(prinds, trids, root) {
  .add <- function(i) {
    nd <- vector("list", length = 4)
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
  tree <- new('TreeMan', ndlst = ndlst, root = root,
              wtxnyms = FALSE, ndmtrx = NULL,
              prinds = prinds, tinds = tinds)
  treeman::checkNdlst(ndlst, root)
  treeman::fastCheckTreeMan
  treeman::updateSlts(tree)
}

#' @name trid_to_txid
#' @title Get taxonomic ID from tree ID
#' @description Return taxonomic ID for taxonomic
#' tree ID in taxonomic dictionary.
#' @return Character or vector
#' @param id trid
#' @param dictionary DictionaryTaxon
#' @noRd
trid_to_txid <- function(trid, dictionary) {
  dictionary@txids[match(trid, dictionary@indx)]
}

#' @name txid_to_trid
#' @title Get tree ID from taxonomic ID
#' @description Return taxonomic tree ID from taxonomic
#' ID in dictionary.
#' @return Character or vector
#' @param id txid
#' @param dictionary DictionaryTaxon
#' @noRd
txid_to_trid <- function(txid, dictionary) {
  as.character(dictionary@indx[match(txid, dictionary@txids)])
}

#' @name txid_singletons_get
#' @title Get equivalent singleton IDs
#' @description For a given taxonomic ID, return all the
#' equivalent singleton IDs.
#' @return Character vector
#' @param txid txid
#' @param dictionary DictionaryTaxon
#' @details phylotaR drops all singleton nodes from
#' the NCBI taxonomy upon initiation. To retrieve all
#' the singleton IDs for an ID use this function.
#' Singletons are defined as nodes in the hierarchy
#' that do not split.
#' @noRd
txid_singletons_get <- function(txid, dictionary) {
  trid <- as.character(dictionary@indx[match(txid, dictionary@txids)])
  sngltns <- dictionary@txids[trid == dictionary@indx]
  sngltns[sngltns != txid]
}

#' @name txid_rank_get
#' @title Get rank
#' @description Look-up taxonomic rank from dictionary.
#' @return Character
#' @param txid txid
#' @param dictionary DictionaryTaxon
#' @noRd
txid_rank_get <- function(txid, dictionary) {
  dictionary@rcrds[[txid]]@rnk
}

#' @name txid_descendants_get
#' @title Get all descendants
#' @description Look-up all taxonomic descendants of a node from
#' taxonomic dictionary.
#' @return Vector
#' @param txid txid
#' @param dictionary DictionaryTaxon
#' @param direct T/F, return only direct descedants?
#' @noRd
txid_descendants_get <- function(txid, dictionary, direct=FALSE) {
  trid <- txid_to_trid(txid = txid, dictionary = dictionary)
  if (direct) {
    tr_ptids <- treeman::getNdSlt(tree = dictionary@txtr,
                                  slt_nm = 'ptid', id = trid)
  } else {
    tr_ptids <- treeman::getNdPtids(tree = dictionary@txtr, id = trid)
  }
  
  trid_to_txid(trid = tr_ptids, dictionary = dictionary)
}

#' @name txid_parent_get
#' @title Get taxonomic parent
#' @description Look-up MRCA of taxonomic id(s) from taxonomic dictionary
#' @return Character
#' @param txid txid(s)
#' @param dictionary DictionaryTaxon
#' @noRd
txid_parent_get <- function(txid, dictionary) {
  trid <- txid_to_trid(txid = txid, dictionary = dictionary)
  if (length(trid) > 1) {
    res <- treeman::getPrnt(tree = dictionary@txtr,
                            ids = unique(trid))
  } else {
    res <- treeman::getNdSlt(tree = dictionary@txtr,
                             slt_nm = 'prid', id = trid)
  }
  trid_to_txid(trid = res, dictionary = dictionary)
}
