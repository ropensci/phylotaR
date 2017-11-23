#' @name genTDObj
#' @title Generate R object from NCBI taxonomy dump
#' @description Generates an interrogable R object
#' from the NCBI taxonomy dump using the \code{CHNOSZ}
#' library. The returned object is a list of taxonomic nodes
#' and names.
#' @param wd Working directory
#' @details Object will be cached.
#' @export
# @Hannes, instead of global env vars, I'm using the cache tools
# to save progress along the way.
genTDObj <- function(wd) {
  if(!chkObj(wd, 'tdobj')) {
    cat('Reading from taxonomy dump ....\n')
    if(!file.exists(file.path(wd, 'NCBI'))) {
      stop('No NCBI folder in `wd`.')
    }
    nds <- CHNOSZ::getnodes(file.path(wd, 'NCBI'))
    nms <- CHNOSZ::getnames(file.path(wd, 'NCBI'))
    tdobj <- list('nds'=nds,
                  'nms'=nms)
    svObj(wd=wd, obj=tdobj, nm='tdobj')
    cat('Done.')
  } else {
    tdobj <- ldObj(wd, 'tdobj')
  }
  tdobj
}

#' @name getMngblIds
#' @title Identify manageable taxonomic node IDs
#' @description Given a root \code{txid}, return a set of
#' taxonomic node IDs that only have up to a maximum number
#' of descendants.
#' @param txid Root taxonomic ID(s), vector or character
#' @param td_nds 'nds' element from the \code{tdobj}
#' @param mx_dscndnts Maximum number of descendants for
#' 'manageable' taxonomic node, numeric.
#' @param tmout Number of seconds before discarding node,
#' numeric.
#' @param verbose print progress to screen? T/F
#' @details Returns a \code{list}.
#' If there are more descendants than mx_dscndnts or
#' retrieveing the number takes longer than \code{tmout}
#' seconds, the taxonomic ID is discarded and its children
#' are added to the queue.
#' @export
# @Hannes: this is functional eqv to get.manageable.node.set
# @Hannes: why do we need the timeout? Surely, mx_dscndnts
# is sufficient?
getMngblIds <- function(txid, td_nds,
                        mx_dscndnts=10000,
                        tmout=10,
                        verbose=FALSE) {
  queue <- txid
  mngbl_ids <- vector()
  ndscndnts <- vector()
  rjctd_ids <- vector()
  tot <- 0
  while(length(queue) > 0) {
    id <- head(queue, 1)
    queue <- tail(queue, length(queue)-1)
    n <- .evlTmLmt(nDscndnts(id, td_nds),
                   cpu=tmout)
    if (is.null(n) || n > mx_dscndnts) {
      queue <- c(queue, children(id, td_nds))
      rjctd_ids <- c(rjctd_ids, id)
      .cp(v=verbose, "Taxon [", id,
          "] has too many descendants or tmout 
reached counting descendants. Processing child taxa.")
    } else {
      mngbl_ids <- c(mngbl_ids, id)
      ndscndnts <- c(ndscndnts, n)
      .cp(v=verbose, "Taxon [", id,
          "] has maneagable number of descendants [",
          n, '].')
      tot <- tot + n
      .cp(v=verbose, "Current number of nodes to be
processed [", tot, "]")
    }
  }
  list('mngbl_ids'=mngbl_ids, 'rjctd_ids'=rjctd_ids,
       'ndscndnts'=ndscndnts)
}

#' @name nDscndnts
#' @title Count descendants
#' @description Count the number of children
#' descending from a node in the NCBI taxonomy
#' dump.
#' @param id Taxonomic ID
#' @param td_nds NCBI taxonomic nodes
#' @export
nDscndnts <- function(id, td_nds) {
  queue <- id
  res <- 0
  while (length(queue) > 0) {
    res <- res + length(queue)
    newqueue <- suppressWarnings(foreach(i=seq_along(queue),
                        .combine=c) %dopar% {
                          getKids(queue[i], td_nds)
                          })
    queue <- newqueue
  }
  return(res-1)
}

#' @name getKids
#' @title Return descendent IDs
#' @description Return vector of descendent IDs
#' from NCBI taxonomic node
#' @export
# @Hannes: eqv of children()
getKids <- function(id, td_nds) {
  td_nds$id[which(td_nds$parent==id)]
}