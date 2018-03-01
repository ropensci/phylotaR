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