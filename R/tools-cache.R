# TODO: lots of duplicated lines, break into small functions.
#' @name progress_init
#' @title Initialise progress list in cache
#' @description Creates a progress list recording each stage run in
#' cache.
#' @param wd Working directory
#' @family run-private
#' @return NULL
progress_init <- function(wd) {
  prgrss <- rep(FALSE, 4)
  names(prgrss) <- c('taxise', 'download', 'cluster', 'cluster2')
  fl <- file.path(wd, 'cache', 'progress.RData')
  saveRDS(object = prgrss, file = fl)
}

#' @name progress_save
#' @title Save current progress
#' @description Stores the pipeline progress in the cache.
#' @param wd Working directory
#' @param stg Stage
#' @family run-private
#' @return NULL
progress_save <- function(wd, stg) {
  fl <- file.path(wd, 'cache', 'progress.RData')
  prgrss <- readRDS(file = fl)
  prgrss[[stg]] <- TRUE
  saveRDS(object = prgrss, file = fl)
}

#' @name progress_reset
#' @title Reset progress
#' @description Reset progress to an earlier completed stage.
#' @details For example, resetting the progress to 'download'
#' mark stages 'download', 'cluster' and 'cluster2' as un-run.
#' @param wd Working directory
#' @param stg Stage to which the pipeline will be reset
#' @family run-private
#' @return NULL
progress_reset <- function(wd, stg) {
  fl <- file.path(wd, 'cache', 'progress.RData')
  prgrss <- readRDS(file = fl)
  i <- which(names(prgrss) == stg)
  prgrss[i:length(prgrss)] <- FALSE
  saveRDS(object = prgrss, file = fl)
}

#' @name progress_read
#' @title Read the progress from cache
#' @description Return the last completed stage using the cache.
#' @param wd Working directory
#' @family run-private
#' @return stage name, character, or NA is complete
progress_read <- function(wd) {
  fl <- file.path(wd, 'cache', 'progress.RData')
  prgrss <- readRDS(file = fl)
  names(prgrss)[sum(prgrss) + 1]
}

#' @name cache_setup
#' @title Set-up a cache
#' @description Creates a cache of parameters in the wd.
#' @template ps
#' @param ovrwrt Overwrite existing cache? Default FALSE.
#' @details Warning: overwriting with this function will delete the
#' existing cache.
#' @family run-private
#' @return NULL
cache_setup <- function(ps, ovrwrt = FALSE) {
  d <- file.path(ps[['wd']], 'cache')
  if (file.exists(d)) {
    if (!ovrwrt) {
      stop('Cache already exists, ovrwrt = FALSE.')
    } else {
      cache_rm(ps[['wd']])
    }
  }
  dir.create(d)
  dir.create(file.path(d, 'ncbi'))
  dir.create(file.path(d, 'ncbi', 'search'))
  dir.create(file.path(d, 'ncbi', 'search', 'taxonomy'))
  dir.create(file.path(d, 'ncbi', 'search', 'nucleotide'))
  dir.create(file.path(d, 'ncbi', 'search', 'nuccore'))
  dir.create(file.path(d, 'ncbi', 'fetch'))
  dir.create(file.path(d, 'ncbi', 'fetch', 'taxonomy'))
  dir.create(file.path(d, 'ncbi', 'fetch', 'nucleotide'))
  dir.create(file.path(d, 'ncbi', 'fetch', 'nuccore'))
  dir.create(file.path(d, 'ncbi', 'summary'))
  dir.create(file.path(d, 'ncbi', 'summary', 'taxonomy'))
  dir.create(file.path(d, 'ncbi', 'summary', 'nucleotide'))
  dir.create(file.path(d, 'ncbi', 'summary', 'nuccore'))
  dir.create(file.path(d, 'clstrs'))
  dir.create(file.path(d, 'sqs'))
  dir.create(file.path(d, 'blast'))
  saveRDS(object = ps, file = file.path(d, 'prmtrs.RData'))
}

#' @name parameters_load
#' @title Load parameters from cache
#' @description Parameters are held in cache, use this function to
#' load parameters set for a wd.
#' @param wd Working directory
#' @family run-private
#' @return Parameters list
parameters_load <- function(wd) {
  fl <- file.path(wd, 'cache', 'prmtrs.RData')
  if (file.exists(fl)) {
    ps <- readRDS(file = fl)
  } else {
    stop('Cache does not exist.')
  }
  ps
}

#' @name obj_check
#' @title Check if an object exists
#' @description Check if an object exists in the cache.
#' @param wd Working directory
#' @param nm Object name
#' @family run-private
#' @return T/F
obj_check <- function(wd, nm) {
  fl <- file.path(wd, 'cache', paste0(nm, '.RData'))
  file.exists(fl)
}

#' @name obj_load
#' @title Load a named object from the cache
#' @description Loads an object from the cache as stored by
#' \code{obj_save}.
#' @param wd Working directory
#' @param nm Object name
#' @family run-private
#' @return object, multiple formats possible
obj_load <- function(wd, nm) {
  fl <- file.path(wd, 'cache', paste0(nm, '.RData'))
  if (!file.exists(fl)) {
    stop(paste0('No [', nm, '] in [',
                wd, '] cache.'))
  }
  readRDS(file = fl)
}

#' @name obj_save
#' @title Save a named object in the cache
#' @description Save an object in the cache that can be loaded by
#' \code{obj_load}.
#' @param wd Working directory
#' @param obj Object
#' @param nm Object name
#' @family run-private
#' @return NULL
obj_save <- function(wd, obj, nm) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    dir.create(d)
  }
  fl <- file.path(wd, 'cache', paste0(nm, '.RData'))
  saveRDS(object = obj, file = fl)
}

#' @name cache_rm
#' @title Delete a cache
#' @description Deletes a cache from a wd.
#' @param wd Working directory
#' @family run-private
#' @return NULL
cache_rm <- function(wd) {
  d <- file.path(wd, 'cache')
  if (file.exists(d)) {
    unlink(d, recursive = TRUE)
  }
}

#' @name sqs_save
#' @title Save sequences to cache
#' @description Saves sequences downloaded
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @param sqs Sequences
#' @details Used within the \code{dwnld} function. Saves
#' sequence data by txid in cache.
#' @family run-private
#' @return NULL
sqs_save <- function(wd, txid, sqs) {
  # TODO: avoid overwriting
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'sqs')
  if (!file.exists(d)) {
    dir.create(d)
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  saveRDS(object = sqs, file = fl)
}

#' @name sqs_load
#' @title Load sequences from cache
#' @description Load sequences downloaded by \code{dwnld} function.
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @family run-private
#' @return SeqArc
sqs_load <- function(wd, txid) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'sqs')
  if (!file.exists(d)) {
    stop('`sqs` not in cache. Have you run the download stage?')
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  if (!file.exists(fl)) {
    stop(paste0('[', txid, '] not in `sqs` of cache.'))
  }
  readRDS(file = fl)
}

#' @name sids_check
#' @title Check if sids exist
#' @description Check if sids are already downloaded for a txid.
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @family run-private
#' @return T/F
sids_check <- function(wd, txid) {
  d <- file.path(wd, 'cache', 'sids')
  fl <- file.path(d, paste0(txid, '.RData'))
  file.exists(fl)
}

#' @name sids_save
#' @title Save sids to cache
#' @description Saves sids downloaded
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @param sids sids
#' @family run-private
#' @return NULL
sids_save <- function(wd, txid, sids) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'sids')
  if (!file.exists(d)) {
    dir.create(d)
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  saveRDS(object = sids, file = fl)
}

#' @name sids_load
#' @title Load sids from cache
#' @description Load sids downloaded by \code{sids_get} function.
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @family run-private
#' @return vector of sids
sids_load <- function(wd, txid) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'sids')
  if (!file.exists(d)) {
    stop('`sids` not in cache. Have you run the download stage?')
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  if (!file.exists(fl)) {
    stop(paste0('[', txid, '] not in `sids` of cache.'))
  }
  readRDS(file = fl)
}

#' @name clstrs_save
#' @title Save clusters to cache
#' @description Saves clusters generated by \code{clstr_all} to cache.
#' @param wd Working directory
#' @param txid Taxonomic ID, numeric
#' @param clstrs cluster list
#' @family run-private
#' @return NULL
clstrs_save <- function(wd, txid, clstrs) {
  d <- file.path(wd, 'cache')
  if (!file.exists(d)) {
    stop('Cache does not exist.')
  }
  d <- file.path(d, 'clstrs')
  if (!file.exists(d)) {
    dir.create(d)
  }
  fl <- file.path(d, paste0(txid, '.RData'))
  saveRDS(object = clstrs, file = fl)
}

#' @name ncbicache_load
#' @title Retrieve cached NCBI query
#' @description Run this function to load cached NCBI queries.
#' @param fnm NCBI Entrez function name
#' @param args Args used for function
#' @param wd Working directory
#' @family run-private
#' @return rentrez result
ncbicache_load <- function(fnm, args, wd) {
  flpth <- file.path(wd, 'cache', 'ncbi',
                     fnm, args[['db']])
  fldctnry_pth <- file.path(flpth, 'fldctnry.RData')
  if (file.exists(fldctnry_pth)) {
    fldctnry <- readRDS(fldctnry_pth)
  } else {
    return(NULL)
  }
  if (fnm == 'search') {
    id <- paste0('TERM_', args[['term']], '_RETSTART_',
                 args[['retstart']], '_RETMAX_', args[['retmax']])
    pull <- vapply(X = fldctnry, FUN = function(x) id == x,
                   FUN.VALUE = logical(1))
  } else {
    id <- args[['id']]
    pull <- vapply(X = fldctnry, FUN.VALUE = logical(1),
                   FUN = function(x) all(id %in% x) & all(x %in% id))
  }
  if (sum(pull) == 0) {
    return(NULL)
  }
  flnm <- paste0(which(pull), '.RData')
  flpth <- file.path(flpth, flnm)
  readRDS(file = flpth)
}

#' @name ncbicache_save
#' @title Save NCBI query result to cache
#' @description Run whenever NCBI queries are made to save results in
#' cache in case the pipeline is run again.
#' @param fnm NCBI Entrez function name
#' @param args Args used for function
#' @param wd Working directory
#' @param obj NCBI query result
#' @family run-private
#' @return NULL
ncbicache_save <- function(fnm, args, wd, obj) {
  flpth <- file.path(wd, 'cache', 'ncbi', fnm, args[['db']])
  fldctnry_pth <- file.path(flpth, 'fldctnry.RData')
  if (file.exists(fldctnry_pth)) {
    fldctnry <- readRDS(fldctnry_pth)
  } else {
    fldctnry <- list()
  }
  if (fnm == 'search') {
    # construct unique identifier from args
    id <- paste0('TERM_', args[['term']], '_RETSTART_',
                 args[['retstart']], '_RETMAX_', args[['retmax']])
  } else {
    # ids are sufficient
    id <- args[['id']]
  }
  fldctnry[[length(fldctnry) + 1]] <- id
  flnm <- paste0(length(fldctnry), '.RData')
  flpth <- file.path(flpth, flnm)
  saveRDS(object = obj, file = flpth)
  saveRDS(object = fldctnry, file = fldctnry_pth)
}

#' @name blastcache_load
#' @title Load BLAST results from cache
#' @description Run to load cached BLAST results.
#' @param sids Sequence IDs
#' @param wd Working dir
#' @family run-private
#' @return blast_res data.frame or NULL
blastcache_load <- function(sids, wd) {
  # fldctnry contains all the IDs of the sequences used in a BLAST
  fldctnry_pth <- file.path(wd, 'cache', 'blast', 'fldctnry.RData')
  if (file.exists(fldctnry_pth)) {
    fldctnry <- readRDS(fldctnry_pth)
  } else {
    return(NULL)
  }
  pull <- vapply(fldctnry, function(x) all(sids %in% x), logical(1))
  if (sum(pull) == 0) {
    return(NULL)
  }
  if (sum(pull) == 0) {
    flnm <- paste0(which(pull), '.RData')
    flpth <- file.path(wd, 'cache', 'blast', flnm)
    return(readRDS(file = flpth))
  }
  pssbls <- which(pull)
  lngs <- vapply(fldctnry[pssbls], length, integer(1))
  flnm <- paste0(pssbls[which.min(lngs)], '.RData')
  flpth <- file.path(wd, 'cache', 'blast', flnm)
  blst_rs <- readRDS(file = flpth)
  pull <- blst_rs[['query.id']] %in% sids &
    blst_rs[['subject.id']] %in% sids
  blst_rs[pull, ]
}

#' @name blastcache_save
#' @title Save BLAST results to cache
#' @description Run whenever local BLAST runs are made to save results
#' in cache in case the pipeline is run again.
#' @param sids Sequence IDs
#' @param wd Working dir
#' @param obj BLAST result
#' @family run-private
#' @return NULL
blastcache_save <- function(sids, wd, obj) {
  # fldctnry contains all the IDs of the sequences used in a BLAST
  fldctnry_pth <- file.path(wd, 'cache', 'blast', 'fldctnry.RData')
  if (file.exists(fldctnry_pth)) {
    fldctnry <- readRDS(fldctnry_pth)
  } else {
    fldctnry <- list()
  }
  fldctnry[[length(fldctnry) + 1]] <- sids
  flnm <- paste0(length(fldctnry), '.RData')
  flpth <- file.path(wd, 'cache', 'blast', flnm)
  saveRDS(object = obj, file = flpth)
  saveRDS(object = fldctnry, file = fldctnry_pth)
}
