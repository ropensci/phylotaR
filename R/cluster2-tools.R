#' @name clstrClstrs
#' @title Cluster sets of clusters identified in cluster stage
#' @description Loads cluster sets from cache. Extracts seed sequences
#' and runs all-v-all BLAST of seeds to identify sister clusters.
#' Sisters are then merged. An object of all sequences and clusters
#' is then saved in cache.
#' @param ps Parameters
#' @export
clstrClstrs <- function(ps) {
  info(lvl=1, ps=ps, 'Loading clusters ...')
  clstrpth <- file.path(ps[['wd']], 'cache', 'clstrs')
  sqpth <- file.path(ps[['wd']], 'cache', 'sqs')
  clstrfls <- list.files(clstrpth)
  sqfls <- list.files(sqpth)
  seeds <- NULL
  all_sqs <- all_clstrs <- list()
  for(i in seq_along(clstrfls)) {
    clstrs <- readRDS(file.path(clstrpth, clstrfls[i]))
    sqs <- readRDS(file.path(sqpth, sqfls[i]))
    all_sqs <- c(all_sqs, sqs)
    all_clstrs <- c(all_clstrs, clstrs)
  }
  if(length(clstrfls) > 1) {
    info(lvl=1, ps=ps, 'Done. Cluster-clustering ...')
    seeds <- getSeedSqs(clstrs=all_clstrs, sqs=all_sqs)
    blst_res <- blstSeeds(sqs=seeds, ps=ps)
    info(lvl=1, ps=ps, 'Done. Merging ...')
    jnd_clstrs <- jnClstrs(blst_res=blst_res, ps=ps,
                           seed_ids=names(seeds),
                           all_clstrs=all_clstrs)
    mrg_clstrs <- mrgClstrs(jnd_clstrs=jnd_clstrs)
    all_clstrs <- c(all_clstrs, mrg_clstrs)
  } else {
    info(lvl=1, ps=ps,
         'Done. Only one cluster set -- skipping cluster^2')
  }
  info(lvl=1, ps=ps, 'Done. Renumbering clusters ...')
  all_clstrs <- rnmbrClstrs(clstrs=all_clstrs)
  info(lvl=1, ps=ps, 'Done. Saving ...')
  svObj(wd=ps[['wd']], obj=list('clstrs'=all_clstrs,
                                'sqs'=all_sqs),
        nm='clstrs_sqs')
}

#' @name blstSeeds
#' @title BLAST seed sequences
#' @description Runs all-v-all blast for seed sequences.
#' @param sqs All seed sequences to be BLASTed
#' @param ps Parameters
#' @export
blstSeeds <- function(sqs, ps) {
  info(lvl=2, ps=ps, "BLASTing all vs all for [",
       length(sqs), "] sequences")
  dbfl <- 'seeds-db.fa'
  outfl <- 'seeds-db-blastout.txt'
  mkBlstDB(sqs, dbfl=dbfl, ps=ps)
  blst_rs <- blstN(dbfl=dbfl, outfl=outfl, ps=ps)
  info(lvl=2, ps=ps, "Done. Number of BLAST results [",
       nrow(blst_rs), "]")
  blst_rs
}

#' @name getSeedSqs
#' @title Extract seed sequences from clusters
#' @description Returns seed sequences of all clusters
#' @param clstrs All clusters
#' @param sqs All seed sequences to be BLASTed
#' @export
getSeedSqs <- function(clstrs, sqs) {
  sids <- vapply(clstrs, function(x) x[['seed_gi']], '')
  sqs[sids]
}

#' @name jnClstrs
#' @title Join clusters for merging
#' @description Uses seed sequence BLAST results and IDs
#' to join clusters identifed as sisters into single clusters.
#' Resulting object is of joined clusters, merging is required
#' to reformat the clusters for subsequent analysis.
#' @param blst_rs Seed sequence BLAST results
#' @param seed_ids Seed sequence IDs
#' @param all_clstrs List of all clusters
#' @param ps Parameters
#' @export
jnClstrs <- function(blst_rs, seed_ids, all_clstrs, ps) {
  join <- function(x) {
    pull <- seed_ids %in% x[['gis']]
    jnd_clstr <- all_clstrs[pull]
    nms <- names(jnd_clstr[[1]])
    clstr <- lapply(nms, function(nm)
      unlist(lapply(jnd_clstr, function(cl) cl[[nm]])))
    names(clstr) <- nms
    clstr[['seed_gi']] <- x[['seed_gi']]
    clstr
  }
  pull <- blst_rs[['query.id']] != blst_rs[['subject.id']] &
    blst_rs[['qcovs']] > ps[['mncvrg']]
  blst_rs <- blst_rs[pull, ]
  clstrs <- clstrBlstRs(blst_rs=blst_rs)
  lapply(clstrs, join)
}

#' @name mrgClstrs
#' @title Merge joined clusters
#' @description Takes a list of joined clusters and computes
#' each data slot to create a single merged cluster.
#' @param jnd_clstrs List of joined clusters
#' @export
mrgClsts <- function(jnd_clstrs) {
  for(i in seq_along(jnd_clstrs)) {
    # TODO: look up parent of both? Might require internet
    jnd_clstrs[[i]][['ti_root']] <- NA
    jnd_clstrs[[i]][['ci']] <- NA
    jnd_clstrs[[i]][['ci_anc']] <- NA
    jnd_clstrs[[i]][['cl_type']] <- 'merged'
    jnd_clstrs[[i]][['n_gi']] <-
      sum(jnd_clstrs[[i]][['n_gi']])
    jnd_clstrs[[i]][['n_ti']] <-
      sum(jnd_clstrs[[i]][['n_ti']])
    jnd_clstrs[[i]][['n_gen']] <-
      sum(jnd_clstrs[[i]][['n_gen']])
    jnd_clstrs[[i]][['n_child']] <-
      sum(jnd_clstrs[[i]][['n_child']])
    jnd_clstrs[[i]][['unique_id']] <-
      paste0(jnd_clstrs[[i]][['unique_id']], collapse='|')
    jnd_clstrs[[i]][['MinLength']] <-
      min(jnd_clstrs[[i]][['MinLength']])
    jnd_clstrs[[i]][['MaxLength']] <-
      max(jnd_clstrs[[i]][['MaxLength']])
  }
  jnd_clstrs
}

#' @name rnmbrClstrs
#' @title Renumbers cluster IDs
#' @description Returns a list of clusters with
#' ID (cis) determined by the number of sequences
#' in each cluster.
#' @param clstrs List of clusters
#' @export
rnmbrClstrs <- function(clstrs) {
  nsqs <- vapply(clstrs, function(x) x[['n_gi']],
                 numeric(1))
  ord <- order(nsqs, decreasing=TRUE)
  clstrs <- clstrs[ord]
  for(i in seq_along(clstrs)) {
    clstrs[[i]][['ci']] <- i - 1
  }
  clstrs
}
