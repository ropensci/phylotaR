#' @name mkBlstDB
#' @title Make a BLAST database
#' @description Generate BLAST database in wd 
#' for given sequences.
#' @param sqs Sequences
#' @param dbfl Outfile for database
#' @param wd Working directory
#' @export
mkBlstDB <- function(sqs, dbfl, wd, verbose=FALSE) {
  if(length(sqs) < 2) {
    stop('Need more than 2 sequences for BLAST.')
  }
  prmtrs <- ldPrmtrs(wd)
  .cp(v=verbose, "Making blast database for [",
      length(sqs), "] sequences")
  blst_d <- file.path(wd, 'blast')
  if(!file.exists(blst_d)) {
    dir.create(blst_d)
  }
  fl <- file.path(blst_d, dbfl)
  file.create(fl)
  for(gi in names(sqs)) {
    write(paste0("> ", gi, "\n", sqs[[as.character(gi)]][['seq']]),
          file=fl, append=TRUE)
  }
  cmd <- paste(prmtrs[['mkblstdb']], '-in', fl, '-dbtype nucl')
  system(cmd) == 0 || stop('Command ', cmd, 'did not return 0')
  # Check success
  extensions <- c('nhr', 'nin', 'nsq')
  fnames <- sapply(extensions, function(e) paste0(fl, '.', e))
  if(!all(sapply(fnames, file.exists))) {
    stop('Command [', cmd, '] did not produce output files [',
         paste(fnames), ']')
  }
  fl
}

#' @name blstN
#' @title BLAST all vs all
#' @description Use \code{blastn} to BLAST all vs all using
#' a BLAST database
#' @param dbfl Database file
#' @param outfl Output file
#' @param wd Working directory
#' @param eval_ctoff Minimum E-value
#' @param verbose Verbose? T/F
#' @export
blstN <- function(dbfl, outfl, wd, eval_ctoff=1.0e-10,
                  verbose=FALSE) {
  prmtrs <- ldPrmtrs(wd)
  blst_d <- file.path(wd, 'blast')
  if(!file.exists(blst_d)) {
    stop('No `blast` dir in wd.')
  }
  outfl <- file.path(blst_d, outfl)
  dbfl <- file.path(blst_d, dbfl)
  if(!file.exists(dbfl)) {
    stop(paste0('[', dbfl, '] does not exist. ',
                'Are you sure you ran `mkBlstDB`?'))
  }
  # TODO: We don't really need all these columns...
  outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'"
  # Disable DUST filtering, limit to same-strang matching
  cmd <- paste(prmtrs[['blstn']], '-query', dbfl, '-db', dbfl, '-outfmt',
               outfmt, '-dust no -strand plus -evalue', eval_ctoff,
               '-out', outfl)
  system(cmd) == 0 || stop(paste0('Command [', cmd, '] did not return 0'))
  if(!file.exists(outfl)) {
    .cp(v=verbose, "No BLAST output, returning NULL")
    return(NULL)
  }
  if(file.info(outfl)[['size']]==0) {
    .cp(v=verbose, "No BLAST output, returning NULL")
    return(NULL)
  }
  blst_rs <- read.table(outfl)
  colnames(blst_rs) <- c('query.id', 'subject.id', 'identity',
                         'alignment.length', 'mismatches',
                         'gap.opens', 'q.start', 'q.end',
                         's.start', 's.end', 'evalue',
                         'bit.score', 'qcovs', 'qcovhsp')
  blst_rs
}

#' @name filter.blast.results
#' @title Filter BLAST results
#' @description TODO
#' @details
#' @export
#' @examples
#' # TODO
filter.blast.results <- function(blast.results, seqs, min.coverage=0.51) {
  # TODO remove the plyr dependency '.'
  blast.dt <- data.table::data.table(blast.results)
  
  ## First filter out HSPs with lower than min.coverage
  ## Get query-subject pairs with too low coverage
  ## Remove both, query-subject and subject-query pair of low coverage hits!!
  ## Otherwise, we will end up with clusters with uneven sequence lengths!
  ## TODO does it make a difference to filter for qcovhsp ??
  remove.dt <- blast.dt[qcovs/100 < min.coverage, .(query.id, subject.id)][, rbind(.SD, rev(.SD),
                                                                                   use.names=FALSE)]
  blast.dt.filtered <- blast.dt[!remove.dt, on=names(remove.dt)]
  cat("Removed", nrow(blast.dt) - nrow(blast.dt.filtered),
      "of the", nrow(blast.dt), "BLAST hits due to insufficient coverage\n")
  
  ## XXX Don't need to collapse when doing single-linkeage clustering!!
  ## collapse HSPs such that we end up with unique query-subject pairs
  ## TODO The xx is a hack.. how do I collapse by two column values and
  ## only get the column back??
  ## result.subset <- setkey(blast.dt.filtered, query.id,
  ## subject.id)[, .(xx=0),.(query.id, subject.id)]
  result <- as.data.frame(blast.dt.filtered[,.(query.id, subject.id)])
  return (result)
}