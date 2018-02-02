#' @name mkBlstDB
#' @title Make a BLAST database
#' @description Generate BLAST database in wd 
#' for given sequences.
#' @param sqs Sequences
#' @param dbfl Outfile for database
#' @export
mkBlstDB <- function(sqs, dbfl, ps) {
  if(length(sqs) < 2) {
    error(ps=ps, 'Need more than 2 sequences for BLAST.')
  }
  info(lvl=3, ps=ps,
       "Making blast database for [", length(sqs), "] sequences")
  blst_d <- file.path(ps[['wd']], 'blast')
  if(!file.exists(blst_d)) {
    dir.create(blst_d)
  }
  fl <- file.path(blst_d, dbfl)
  file.create(fl)
  for(gi in names(sqs)) {
    write(paste0("> ", gi, "\n", sqs[[as.character(gi)]][['seq']]),
          file=fl, append=TRUE)
  }
  args <- c('-in', fl, '-dbtype nucl')
  res <- .system(command=ps[['mkblstdb']], args=args,
                 stdout=FALSE, stderr=FALSE)
  if(res != 0) {
    error(paste0('Command did not return 0'))
  }
  # Check success
  extensions <- c('nhr', 'nin', 'nsq')
  fnames <- sapply(extensions, function(e) paste0(fl, '.', e))
  if(!all(sapply(fnames, file.exists))) {
    error(ps=ps, 'Command [', cmd, '] did not produce output files [',
          paste(fnames), ']')
  }
  NULL
}

#' @name blstN
#' @title BLAST all vs all
#' @description Use \code{blastn} to BLAST all vs all using
#' a BLAST database
#' @param dbfl Database file
#' @param outfl Output file
#' @export
blstN <- function(dbfl, outfl, ps) {
  blst_d <- file.path(ps[['wd']], 'blast')
  if(!file.exists(blst_d)) {
    error(ps=ps, 'No `blast` dir in wd.')
  }
  outfl <- file.path(blst_d, outfl)
  dbfl <- file.path(blst_d, dbfl)
  if(!file.exists(dbfl)) {
    error(ps=ps, paste0('[', dbfl, '] does not exist. ',
                        'Are you sure you ran `mkBlstDB`?'))
  }
  # TODO: We don't really need all these columns...
  outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'"
  # Disable DUST filtering, limit to same-strang matching
  args <- c('-query', dbfl, '-db', dbfl, '-outfmt',
            outfmt, '-dust no -strand plus -evalue', ps[['mxeval']],
            '-out', outfl)
  res <- .system(command=ps[['blstn']], args=args, stdout=FALSE,
                 stderr=FALSE)
  if(res != 0) {
    error(paste0('Command [', cmd, '] did not return 0'))
  }
  if(!file.exists(outfl)) {
    info(lvl=3, ps=ps, "No BLAST output, returning NULL")
    return(NULL)
  }
  if(file.info(outfl)[['size']]==0) {
    info(lvl=3, ps=ps, "No BLAST output, returning NULL")
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

# FILTER NOTES
# First filter out HSPs with lower than min.coverage
# Get query-subject pairs with too low coverage
# Remove both, query-subject and subject-query pair of low coverage hits!!
# Otherwise, we will end up with clusters with uneven sequence lengths!
# TODO does it make a difference to filter for qcovhsp ??
# XXX Don't need to collapse when doing single-linkeage clustering!!
# collapse HSPs such that we end up with unique query-subject pairs
# TODO The xx is a hack.. how do I collapse by two column values and
# only get the column back??
# result.subset <- setkey(blast.dt.filtered, query.id,
# subject.id)[, .(xx=0),.(query.id, subject.id)]
#' @name fltrBlstRs
#' @title Filter BLAST results
#' @description Can you fill this in Hannes?
#' @param blst_rs BLAST results
#' @export
# @Hannes: do you have any blast results to test this?
fltrBlstRs <- function(blst_rs, ps) {
  pull <- blst_rs[['qcovs']] < ps[['mncvrg']]
  if(any(pull)) {
    # drop all with < mncvrg
    # ensure both qry-sbj and sbj-qry are dropped
    qsids <- blst_rs[pull, c('query.id', 'subject.id')]
    pull <- (blst_rs[['query.id']] %in% qsids[['subject.id']] &
      blst_rs[['subject.id']] %in% qsids[['query.id']]) | pull
  }
  ndrp <- sum(pull)
  info(lvl=3, ps=ps, "Removed [", ndrp, "/",
       nrow(blst_rs), "] BLAST hits due to insufficient coverage")
  if(sum(pull) > 0) {
    blst_rs <- blst_rs[!pull, ]
  }
  blst_rs
}
