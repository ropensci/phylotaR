#' @name calcClstrs
#' @title Calculate clusters for all sequences in WD
#' @description TODO
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
calcClstrs <- function(phylt_nds, ps) {
  # load sequences
  sq_fls <- list.files(file.path(ps[['wd']], 'cache', 'sqs'))
  #foreach(i=seq_along(sq_fls)) %dopar% {
  for(i in seq_along(sq_fls)) {
    sq_fl <- sq_fls[i]
    # TODO: use the cache tool
    sqs <- readRDS(file=file.path(file.path(ps[['wd']], 'cache',
                                            'sqs', sq_fl)))
    txid <- as.numeric(sub('\\.RData', '', sq_fl))
    info(lvl=1, ps=ps, "Working on [id ", txid, "]")
    clstrs <- clstrAll(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
                       ps=ps)
    svClstrs(wd=ps[['wd']], txid=txid, clstrs=clstrs)
    info(lvl=1, ps=ps, "Finished [id ", txid, "] # [",
        i, "/", length(sq_fls), "]")
  }
}

#' @name writeClstr
#' @title Write out PhyLoTa cluster
#' @description Takes PhyLoTa data.frame and writes
#' as .tsv
#' @param phylt_nds PhyLoTa nodes data.frame
#' @details PhyLoTa data.frame must be informed by
#' clustering functions before writing out.
#' @export
writeClstr <- function(txid, cldf, cigidf, sqs, ps) {
  dbpth <- file.path(ps[['wd']], 'dbfiles')
  if(!file.exists(dbpth)) {
    dir.create(dbpth)
  }
  clstrs_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                       '-clusters.tsv'))
  sqs_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                    '-seqs.tsv'))
  ci_gi_fl <- file.path(dbpth, paste0('dbfiles-', txid,
                                      '-ci_gi.tsv'))
  sqdf <- do.call(rbind, lapply(sqs, as.data.frame))
  info(lvl=1, ps=ps, "Taxid", ps[['txid']], ": writing",
       nrow(cldf), "clusters,", nrow(sqdf), "sequences,",
       nrow(cigidf), "ci_gi entries to file")
  write.table(cldf, file=clstrs_fl, append=file.exists(clstrs_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(clstrs_fl))
  write.table(sqdf, file=sqs_fl, append=file.exists(sqs_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(sqs_fl))
  write.table(cigidf, file=ci_gi_fl, append=file.exists(ci_gi_fl),
              quote=FALSE, sep="\t", row.names=FALSE,
              col.names=!file.exists(ci_gi_fl))
}

#' @name clstrAll
#' @title Hierarchically cluster all sequences of a txid
#' @description Identifies all direct and subtree clusters
#' for a taxonomic ID.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of all downloaded sequences
#' @param phylt_nds PhyLoTa table
#' @param ps Parameters
#' @param lvl Log level
#' @export
clstrAll <- function(txid, sqs, phylt_nds, ps, lvl=0) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  all_clstrs <- clstrSbtr(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
                          ps=ps, dds=dds, lvl=lvl+1)
  for(dd in dds) {
    info(lvl=lvl+2, ps=ps, "Processing child [id ", txid,
         "] [", dd, "]")
    dd_clstrs <- clstrAll(txid=dd, phylt_nds=phylt_nds, 
                          sqs=sqs, ps=ps, lvl=lvl+1)
    dd_clstrs <- lapply(dd_clstrs, function(cc) {
      idx <- max(which(sapply(dd_clstrs, function(c) any(
        cc[['gis']] %in% c[['gis']]))))
      cc[['ci_anc']] <- dd_clstrs[[idx]][['ci']]
      dd_clstrs[[idx]][['n_child']] <-
        dd_clstrs[[idx]][['n_child']] + 1
      cc
    })
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
#' @export
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
  info(lvl=lvl+2, ps=ps, "[", length(sids), " sqs]")
  if(length(sids) < 3) {
    info(lvl=lvl+3, ps=ps,
         "Too few sequences, cannot make clusters")
    return(all_clstrs)
  }
  sqs_prt <- sqs[sids]
  sbtr_clstrs <- clstrSqs(txid=txid, sqs=sqs_prt,
                          phylt_nds=phylt_nds,
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
#' @export
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
           phylt_nds=phylt_nds, ps=ps, lvl=lvl)
}

#' @name clstrSqs
#' @title Identify clusters from sequences
#' @description Given a sequence object, this function will generate
#' PhyLoTa-like tables using BLAST.
#' @param txid Taxonomic ID
#' @param sqs Sequence object of sequences to be BLASTed
#' @param phylt_nds PhyLoTa table
#' @param ps Parameters
#' @param typ Direct or Subtree?
#' @export
clstrSqs <- function(txid, sqs, phylt_nds, ps,
                     lvl, typ=c('direct', 'subtree')) {
  typ <- match.arg(typ)
  info(lvl=lvl+1, ps=ps, "BLASTing [", length(sqs@ids),
       " sqs] ....")
  blst_rs <- blstSqs(txid=txid, typ=typ, sqs=sqs, ps=ps, lvl=lvl)
  if(is.null(blst_rs)) {
    return(NULL)
  }
  raw_clstrs <- clstrBlstRs(blst_rs=blst_rs)
  clstrs <- addClstrInf(clstrs=raw_clstrs, phylt_nds=phylt_nds,
                        txid=txid, sqs=sqs, typ=typ)
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
#' @export
# TODO: create separate taxonomy look-up tools
getADs <- function(txid, phylt_nds) {
  dds <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  res <- dds
  for(dd in dds) {
    res <- c(res, getADs(txid=dd, phylt_nds=phylt_nds))
  }
  res
}

#' @name getGnsFrmPhyltNds
#' @title Get genus for txid using PhyLoTa nodes
#' @description Returns genus for a taxonomic node ID
#' @param txid Taxonomic node ID, numeric
#' @param phylt_nds PhyLoTa nodes data.frame
#' @export
getGnsFrmPhyltNds <- function(txid, phylt_nds) {
  # @hannes can you check that this is doing what we want?
  # this always be 1 (below genus) or none (above genus)
  res <- unique(phylt_nds[phylt_nds[['ti_anc']]==txid, 'ti_genus'])
  if(length(res) > 1 || length(res) == 0) {
    res <- NA
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
#' @export
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
#' @export
clstrBlstRs <- function(blst_rs) {
  g <- igraph::graph.data.frame(blst_rs[ ,c("query.id",
                                            "subject.id")],
                                directed=FALSE)
  clstrs <- igraph::clusters(g)
  # filter for phylogenetically informative clusters
  clstrs <- clstrs[['membership']]
  # we will return a list, one entry with sequence IDs
  #  for each cluster
  clstr_lst <- lapply(unique(clstrs), function(x) {
    list(gis=sort(names(clstrs)[which(clstrs==x)]))
  })
  # Get the seed gi, we will chose it to be the sequence
  # in the cluster that has
  # the most hits with the other members in the cluster;
  # i.e. the most connected
  # node in the graph
  degrees <- igraph::degree(g)
  # get seed gis and as field to clusters
  clstr_lst <- lapply(clstr_lst, function(cl){
    idx <- order(degrees[cl[['gis']]], decreasing=TRUE)[1]
    # index of most connected component
    cl[['seed_gi']] <- cl[['gis']][idx]
    cl
  })
  clstr_lst
}

#' @name addClstrInf
#' @title Add taxonomic and sequence info to clusters
#' @description input: List of clstrs
#' @param clstrs Clusters
#' @param phylt_nds PhyLoTa data.frame
#' @param txid Taxonomic node ID
#' @param sqs Sequneces
#' @param typ Subtree of direct?
#' @export
addClstrInf <- function(clstrs, phylt_nds, txid, sqs, typ) {
  # we will take the column names in the database as names in the list
  for(i in 1:length(clstrs)) {
    cl <- clstrs[[i]]
    cl_sqs <- lapply(cl[['gis']], function(gi)sqs[[as.character(gi)]])
    cl[['ti_root']] <- txid
    cl[['ci']] <- i-1
    cl[['cl_type']] <- typ
    cl[['n_gi']] <- length(cl[['gis']])
    # all taxon ids for gis in cluster
    cl[['tis']] <- sapply(cl_sqs, function(x) x@txid)
    cl[['n_ti']] <- length(unique(cl[['tis']]))
    # sequence lengths
    l <- sapply(cl_sqs, function(x) x@nncltds)
    cl[['MinLength']] <- min(l)
    cl[['MaxLength']] <- max(l)
    # get n. genera
    cl[['n_gen']] <- length(unique(sapply(cl[['tis']],
                                          getGnsFrmPhyltNds, phylt_nds)))
    # n_child and ci_anc will be set (or incremented) later,
    # when multiple hierarchies are calculated
    cl[['n_child']] <- 0
    cl[['ci_anc']] <- NA
    # make a unique cluster id consisting of seed gi, taxon id,
    # cluster id, cluster type
    unique_id <- paste0(cl[['seed_gi']], '-', txid, '-',
                        cl[['ci']], '-', cl[['cl_type']])
    cl[['unique_id']] <- unique_id
    clstrs[[i]] <- cl
  }
  clstrs
}

#' @name clstrPhylt
#' @title Reformat cluster obj for PhyLoTa table
#' @description Returns a data frame with columns matchine the fields in
#' PhyLoTA's cluster table
#' @param clstrs Clusters
#' @export
clstrPhylt <- function(clstrs) {
  # remove gi and id fields which won't be part of cluster table
  l <- lapply(clstrs, function(cl)cl[which(!names(cl)%in%c(
    'gis', 'tis', 'unique_id'))])
  # create data frame
  do.call(rbind, lapply(l, as.data.frame))
}

#' @name clstrCiGi
#' @title Reformat cluster obj into CI GI
#' @description TODO
#' @param clstrs Clusters
#' @export
clstrCiGi <- function(clstrs) {
  ## select relevant columns from cluster objects and make data frame    
  do.call(rbind, lapply(clstrs, function(c) {
    data.frame(ti=rep(c[['ti_root']], c[['n_gi']]),
               clustid=rep(c[['ci']], c[['n_gi']]),
               cl_type=rep(c[['cl_type']], c[['n_gi']]),
               gi=c[['gis']],
               ti_of_gi=c[['tis']])
  }))
}
