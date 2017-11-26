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

#' @name genPhylotaNds
#' @title Generate taxonomic nodes in PhyLoTa format
#' @description
#' @param nids Taxonomic node IDs, character vector
#' @param verbose Verbose? T/F
#' @details
#' @export
# @Hannes: this is functional eqv to bulk of nodes.create
genPhylotaNds <- function(nids, verbose=FALSE) {
  phylt_nds <- data.frame('ti'=NA,
                          'ti_anc'=NA,
                          'rank'=NA,
                          'n_gi_node'=NA,
                          'n_gi_sub_nonmodel'=NA,
                          'n_gi_sub_model'=NA,
                          'n_sp_desc'=NA,
                          'n_sp_model'=NA,
                          'n_leaf_desc'=NA,
                          'n_otu_desc'=NA,
                          'ti_genus'=NA,
                          'n_genera'=NA)
  .cp(v=verbose, "Adding [", length(nids),
      "] nodes recursively, in parallel")
  for(i in seq_along(nids)) {
    tax <- nids[i]
    cat("Recursively processing txid ", tax, " # ", i, " / ",
        length(nids), "\n")
    n <- .add.stats(txid=tax, nodes.written=start, 
                    nodes=ncbi.nodes, recursive=TRUE, prmtrs=prmtrs)
    ids.to.write <- n$ti[which(!n$ti %in% start$ti)]
    stopifnot(tax %in% ids.to.write)
    ## only return rows which were not written to file before
    subset <-n[match(ids.to.write, n$ti), ]
    .write.row(subset, ncbi.names, file.name)
    
    cat("Finished processing txid ", tax, " # ", i, " / ",
        length(nids), "\n")
  }
  
  # if length(nids) == 0, file is not initiated
  if(length(nids) == 0) {
    start <- .init.nodes.df()
    write.table(x=start, file=file.name)
  }
  
  ## Add the top nodes non-recursively. This has to be in reversed
  # order to make sure a parent is not
  ## inserted before it's children, since we need the info from
  # the children. Therefore, this must not
  ## happen in parallel!
  cat("Adding ", length(ids.not.processed),
      " top-level nodes non-recursively, sequentially \n")
  for (i in seq_along(ids.not.processed)) {
    tax <- rev(ids.not.processed)[i]
    ##foreach(i=seq_along(ids.not.processed)) %dopar% {
    ## reload nodes that have been written before
    nodes.written <- read.table(file.name, header=TRUE)
    nodes.written <- nodes.written[,c('ti','ti_anc','rank',
                                      'n_gi_node','n_gi_sub_nonmodel',
                                      'n_gi_sub_model',
                                      'n_sp_desc', 'n_sp_model',
                                      'n_leaf_desc', 'n_otu_desc',
                                      'ti_genus',
                                      'n_genera')]
    
    cat("Processing txid ", tax, " # ", i, " / ",
        length(ids.not.processed), "\n")
    n <- .add.stats(tax, nodes.written, ncbi.nodes, recursive=FALSE)
    .write.row(n[match(tax, n$ti),], ncbi.names, file.name)
    cat("Finished processing txid ", tax, " # ",
        i, " / ", length(ids.not.processed), "\n")
  }
}

#' @name getStats
#' @title Generate PhyLoTa statistics for txid
#' @description 
#' @param txid Root taxonomic ID(s), vector or character
#' @param phylt_nds PhyLoTa data.frame
#' @param td_nds 'nds' element from the \code{tdobj}
#' @param td_nms 'nms' element from the \code{tdobj}
#' @param mx_sq_lngth Maximum number of associated sequences, numeric
#' @param mdl_threshold Maximum threshold, numeric.
#' @param verbose Print progress to screen? T/F
#' @param recursive Add stats for kid nodes recursively? T/F
#' @details Returns an updated \code{phylt_nds}.
#' @export
# @Hannes: this is functional eqv to get.stats
getStats <- function(txid, phylt_nds,
                     td_nds, td_nms,
                     mx_sq_lngth,
                     mdl_thrshld,
                     verbose=FALSE,
                     recursive=FALSE) {
  .cp(v=verbose, "Adding species counts for txid [",
      txid, "]")
  rank <- CHNOSZ::getrank(txid, nodes=td_nds)
  parent <- CHNOSZ::parent(txid, nodes=td_nds)
  genus <- getGenus(txid, td_nds)
  stats <- data.frame(ti=txid,
                      ti_anc=parent,
                      rank=rank,
                      n_gi_node=0,
                      n_gi_sub_model=0,
                      n_gi_sub_nonmodel=0,
                      n_leaf_desc=0, # counts itself, for a leaf it is 1
                      n_sp_desc=0, # counts itself, for a species it is 1
                      n_sp_model=0,
                      ti_genus=genus,
                      n_genera=0,
                      n_otu_desc=1)
  if(rank == 'species') {
    stats['n_sp_desc'] <- 1
  }
  n_drctsqs <- nSqs(txid, direct=TRUE,
                    mx_len=mx_sq_lngth,
                    verbose=verbose)
  stats['n_gi_node'] <- n_drctsqs
  kids <- getKids(txid, td_nds=td_nds)
  if(length(kids) == 0) {
    stats['n_leaf_desc'] <- stats['n_leaf_desc'] + 1
  }
  if(n_drctsqs > mdl_thrshld) {
    stats['n_gi_sub_model'] <- n_drctsqs
    stats['n_sp_model'] <- stats['n_sp_model'] + 1
  } else {
    stats['n_gi_sub_nonmodel'] <- n_drctsqs
  }
  for(kid in kids) {
    if(recursive) {
      phylt_nds <- getStats(txid=kid, phylt_nds=phylt_nds,
                            mx_sq_lngth=mx_sq_lngth,
                            mdl_thrshld=mdl_thrshld,
                            td_nds=td_nds, td_nms=td_nms,
                            verbose=verbose, recursive=recursive)
    }
    # TODO: @Hannes, this section worries me. What if a user
    #  supplies phylt_nds without rows for any kids?
    #  Function crashes.
    if(CHNOSZ::getrank(kid, nodes=td_nds) == 'genus') {
      stats['n_genera'] <- stats['n_genera'] + 1
    }
    cols <- c('n_leaf_desc', 'n_otu_desc', 'n_sp_desc',
              'n_sp_model', 'n_genera')
    kid_nd <- phylt_nds[which(phylt_nds[,'ti']==kid),]
    stats['n_gi_sub_model'] <- stats['n_gi_sub_model'] +
      kid_nd['n_gi_sub_model']
    stats['n_gi_sub_nonmodel'] <- stats['n_gi_sub_nonmodel'] +
      kid_nd["n_gi_sub_nonmodel"]
    stats[cols] <- stats[cols] + phylt_nds[
      which(phylt_nds$ti==kid), cols]
    stopifnot(dim(stats)[1]==1)
  }
  .cp(v=verbose, "Added stats for txid [", txid, "]")
  rbind(phylt_nds, stats[names(phylt_nds)])
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

# TODO: rename -- 'kids' implies tips.
#' @name getKids
#' @title Return descendent IDs
#' @description Return vector of descendent IDs
#' from NCBI taxonomic node
#' @export
# @Hannes: eqv of children()
# TODO: is this really 'children' and not just sub nodes?
getKids <- function(id, td_nds) {
  td_nds$id[which(td_nds$parent==id)]
}

#' @name getGenus
#' @title Return genus ID
#' @description Return vector of descendent IDs
#' from NCBI taxonomic node
#' @export
# @Hannes: eqv of .get.genus()
getGenus <- function(txid, td_nds) {
  lwr_rnks <- c('genus', 'subgenus', 'species group',
                'species subgroup', 'species',
                'subspecies','varietas', 'forma')
  hghr_rnks <- c('superkingdom', 'kingdom', 'subkingdom',
                 'superphylum', 'phylum',
                 'subphylum', 'superclass', 'class',
                 'subclass', 'infraclass', 'superorder',
                 'order', 'suborder', 'infraorder',
                 'parvorder', 'superfamily', 'family',
                 'subfamily', 'tribe', 'subtribe')
  # if rank is already higher than genus,
  #  do not search for it
  crrnt <- CHNOSZ::getrank(txid, nodes=td_nds)
  if(!crrnt  %in% lwr_rnks) {
    return(NA)
  }
  while(crrnt != 'genus') {
    txid <- CHNOSZ::parent(txid, nodes=td_nds)
    crrnt <- CHNOSZ::getrank(txid, nodes=td_nds)
    # some hybrids between genera can have a
    #  subfamily as a parent...
    if(crrnt %in% hghr_rnks) {
      break
    }
  }
  return(txid)
}
