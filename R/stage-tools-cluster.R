#' @name clstrAll
#' @title Hierarchically cluster all sequences of a txid
#' @description Identifies all direct and subtree clusters
#' for a taxonomic ID.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param phylt_nds PhyLoTa table
#' @param ps Parameters
#' @param lvl Log level
clstrAll <- function(txid, sqs, phylt_nds, ps, lvl=0) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  all_clstrs <- clstrSbtr(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
                          ps=ps, dds=dds, lvl=lvl+1)
  for(dd in dds) {
    info(lvl=lvl+2, ps=ps, "Processing [id ", txid,
         "] child [id ", dd, "]")
    dd_clstrs <- clstrAll(txid=dd, phylt_nds=phylt_nds, 
                          sqs=sqs, ps=ps, lvl=lvl+1)
    all_clstrs <- c(all_clstrs, dd_clstrs)
  }
  all_clstrs
}

#' @name clstrSbtr
#' @title Cluster all sequences descending from a txid
#' @description Identifies clusters from sequences associated
#' with a txid and all its descendants. Clusters returned by
#' this function will thus be of cl_type 'subtree'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param phylt_nds PhyLoTa table
#' @param dds Vector of direct descendants
#' @param ps Parameters
clstrSbtr <- function(txid, sqs, phylt_nds, dds, ps, lvl) {
  all_clstrs <- list()
  rnks <- as.character(phylt_nds[['rank']])
  rnk <- as.character(rnks[match(txid, phylt_nds[['ti']])])
  info(lvl=lvl+1, ps=ps, "Generating subtree clusters for [id ",
       txid, "(", rnk, ")]")
  if(length(dds) > 0) {
    drct_clstrs <- clstrDrct(txid, ps=ps, phylt_nds=phylt_nds,
                             sqs=sqs, lvl=lvl)
    all_clstrs <- c(all_clstrs, drct_clstrs)
  }
  txids <- getADs(txid=txid, phylt_nds=phylt_nds)
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[which(all_sq_txids %in% as.character(txids))]
  if(length(sids) < 3) {
    info(lvl=lvl+3, ps=ps, "[", length(sids), " sqs]",
         " -- too few sequences, cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sids]
  sbtr_clstrs <- clstrSqs(txid=txid, sqs=sqs_prt,
                          typ='subtree',
                          ps=ps, lvl=lvl)
  c(all_clstrs, sbtr_clstrs)
}

#' @name clstrDrct
#' @title Cluster sequences directly associated with txid
#' @description In GenBank certain sequences may only be associated
#' with a higher level taxon (e.g. genus, family ...). This function
#' generates clusters from these sequences, alone. This function
#' identifies such sequences in the sequence object and generates
#' a list of clusters of cl_type 'Node'.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param phylt_nds PhyLoTa table
#' @param ps Parameters
clstrDrct <- function(txid, sqs, phylt_nds, ps, lvl) {
  all_clstrs <- list()
  rnks <- as.character(phylt_nds[['rank']])
  rnk <- as.character(rnks[match(txid, phylt_nds[['ti']])])
  info(lvl=lvl+1, ps=ps, "Generating direct clusters for [id ",
       txid, "(", rnk, ")]")
  all_sq_txids <- sqs@txids
  sids <- sqs@ids[all_sq_txids %in% as.character(txid)]
  info(lvl=lvl+2, ps=ps, "[", length(sids), " sqs]")
  if(length(sids) < 3) {
    info(lvl=lvl+3, ps=ps,
         "Too few sequences, cannot make clusters")
    return(list())
  }
  sqs_prt <- sqs[sids]
  clstrSqs(txid=txid, sqs=sqs_prt, typ='direct',
           ps=ps, lvl=lvl)
}

#' @name clstrSqs
#' @title Identify clusters from sequences
#' @description Given a sequence object, this function will generate
#' a list of Clstr objects using BLAST
#' @param txid Taxonomic ID
#' @param sqs Sequence object of sequences to be BLASTed
#' @param ps Parameters
#' @param typ Direct or Subtree?
clstrSqs <- function(txid, sqs, ps, lvl,
                     typ=c('direct', 'subtree')) {
  typ <- match.arg(typ)
  info(lvl=lvl+1, ps=ps, "BLASTing [", length(sqs@ids),
       " sqs] ....")
  blst_rs <- blstSqs(txid=txid, typ=typ, sqs=sqs, ps=ps, lvl=lvl)
  if(is.null(blst_rs)) {
    return(NULL)
  }
  clstr_lst <- clstrBlstRs(blst_rs=blst_rs)
  clstrs <- genClstr(clstr_lst=clstr_lst, txid=txid,
                     sqs=sqs, typ=typ)
  info(lvl=lvl+1, ps=ps, "Identified [", length(clstrs),
       "] clusters")
  clstrs
}

#' @name getADs
#' @title Get all descendants
#' @description Return all the taxonomic node IDs descending
#' from given taxonomic ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
# TODO: create separate taxonomy look-up tools
getADs <- function(txid, phylt_nds) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  res <- dds
  for(dd in dds) {
    res <- c(res, getADs(txid=dd, phylt_nds=phylt_nds))
  }
  res
}

#' @name blstSqs
#' @title BLAST All vs All
#' @description Return blast results from
#' BLASTing all vs all for given sequences
#' @param txid Taxonomic node ID, numeric
#' @param typ Cluster type, 'direct' or 'subtree'
#' @param sqs Sequences
#' @param wd Working directory
#' @param verbose Verbose? T/F
blstSqs <- function(txid, typ, sqs, ps, lvl) {
  blst_rs <- ldBlstCch(sqs@ids, wd=ps[['wd']])
  if(is.null(blst_rs)) {
    dbfl <- paste0('taxon-', txid, '-typ-', typ,
                   '-db.fa')
    outfl <- paste0('taxon-', txid, '-typ-', typ,
                    '-blastout.txt')
    mkBlstDB(sqs=sqs, dbfl=dbfl, ps=ps)
    blst_rs <- blstN(dbfl=dbfl, outfl=outfl, ps=ps)
    if(is.null(blst_rs)) {
      blst_rs <- NA
    }
    svBlstCch(sqs@ids, wd=ps[['wd']], obj=blst_rs)
  }
  # TODO: Not so elegant
  if(any(is.na(blst_rs))) {
    return(NULL)
  }
  fltrBlstRs(blst_rs=blst_rs, ps=ps)
}

#' @name clstrBlstRs
#' @title Cluster BLAST Results
#' @description TODO
#' @param blst_rs BLAST results
clstrBlstRs <- function(blst_rs) {
  g <- igraph::graph.data.frame(blst_rs[ ,c("query.id",
                                            "subject.id")],
                                directed=FALSE)
  clstrs <- igraph::clusters(g)
  clstrs <- clstrs[['membership']]
  clstr_lst <- lapply(unique(clstrs), function(x) {
    list('sids'=sort(names(clstrs)[which(clstrs==x)]))
  })
  degrees <- igraph::degree(g)
  clstr_lst <- lapply(clstr_lst, function(cl){
    idx <- order(degrees[cl[['sids']]], decreasing=TRUE)[1]
    # index of most connected component
    cl[['seed']] <- cl[['sids']][idx]
    cl
  })
  clstr_lst
}

#' @name genClstr
#' @title Generate list of Clstrs
#' @description Takes a list of lists of cluster descriptions,
#' returns a list of Clstr objects.
#' @param clstr_lst List of list of cluster descriptions
#' @param txid Taxonomic node ID
#' @param sqs Sequnece records
#' @param typ Subtree of direct?
genClstr <- function(clstr_lst, txid, sqs, typ) {
  clstrs <- vector('list', length=length(clstr_lst))
  for(i in seq_along(clstr_lst)) {
    cl <- clstr_lst[[i]]
    cl_sqs <- sqs[cl[['sids']]]
    clstr <- new('Clstr', sids=cl[['sids']], txids=cl_sqs@txids,
                 typ=typ, prnt=as.character(txid), seed=cl[['seed']])
    clstrs[[i]] <- clstr
  }
  clstrs
}