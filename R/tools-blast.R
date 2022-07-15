#' @name blastdb_gen
#' @title Generate a BLAST database
#' @description Generate BLAST database in wd for given sequences.
#' @param sqs Sequences
#' @param dbfl Outfile for database
#' @template ps
#' @return NULL
#' @family run-private
blastdb_gen <- function(sqs, dbfl, ps) {
  if (length(sqs@sqs) < 2) {
    error(ps = ps, "Need more than 2 sequences for BLAST.")
  }
  blast_d <- file.path(ps[["wd"]], "blast")
  if (!file.exists(blast_d)) {
    dir.create(blast_d)
  }
  fl <- file.path(blast_d, dbfl)
  file.create(fl)
  for (i in seq_along(sqs@ids)) {
    sq <- sqs@sqs[[i]]
    write(paste0("> ", sq@id, "\n", rawToChar(sq@sq)),
      file = fl, append = TRUE
    )
  }
  args <- c("-in", fl, "-dbtype", "nucl")
  info(lvl = 3, ps = ps, "Running makeblastdb")
  res <- cmdln(
    cmd = ps[["mkblstdb"]], args = args,
    lgfl = fl, ps = ps
  )
  if (res != 0) {
    error(ps = ps, paste0("makeblastdb failed to run. Check BLAST log files."))
  }
  # Check success
  extensions <- c("nhr", "nin", "nsq")
  fnames <- vapply(extensions, function(e) paste0(fl, ".", e), character(1))
  if (!all(vapply(fnames, file.exists, logical(1)))) {
    error(ps = ps, "Command did not produce output files [", paste(fnames), "]")
  }
  NULL
}

#' @name blastn_run
#' @title Launch blastn
#' @description Use \code{blastn} to BLAST all-vs-all using a BLAST
#' database.
#' @param dbfl Database file
#' @param outfl Output file
#' @template ps
#' @family run-private
#' @return NULL
blastn_run <- function(dbfl, outfl, ps) {
  blst_d <- file.path(ps[["wd"]], "blast")
  if (!file.exists(blst_d)) {
    error(ps = ps, "No `blast` dir in wd.")
  }
  outfl <- file.path(blst_d, outfl)
  dbfl <- file.path(blst_d, dbfl)
  if (!file.exists(dbfl)) {
    error(ps = ps, paste0("[", dbfl, "] does not exist. "))
  }
  args <- c(
    "-query", dbfl, "-db", dbfl, "-outfmt", ps[["outfmt"]], "-dust",
    "no", "-strand", "plus", "-evalue", ps[["mxevl"]], "-out", outfl
  )
  if (ps[["ncps"]] > 1) {
    args <- c(args, "-num_threads", ps[["ncps"]])
  }
  info(lvl = 3, ps = ps, "Running blastn")
  res <- cmdln(cmd = ps[["blstn"]], args = args, lgfl = dbfl, ps = ps)
  if (res != 0) {
    error(ps = ps, "blastn failed to run. Check BLAST log files.")
  }
  if (!file.exists(outfl)) {
    info(lvl = 3, ps = ps, "No BLAST output, returning NULL")
    return(NULL)
  }
  if (file.info(outfl)[["size"]] == 0) {
    info(lvl = 3, ps = ps, "No BLAST output, returning NULL")
    return(NULL)
  }
  blast_res <- read.table(file = outfl, stringsAsFactors = FALSE)
  colnames(blast_res) <- c(
    "query.id", "subject.id", "identity",
    "alignment.length", "evalue",
    "qcovs", "qcovhsp"
  )
  blast_res
}

#' @name outfmt_get
#' @title Determine 'outformat' format
#' @description Depending on operating system, BLAST may or may not require ""
#' around \code{outfmt}. This function will run a micro BLAST analysis
#' to test. It will return the \code{outfmt} for use in \code{blastn_run}.
#' @template ps
#' @family run-private
#' @return character
#' @export
# Patch for issue 39: https://github.com/ropensci/phylotaR/issues/39
outfmt_get <- function(ps) {
  dbfl <- "testdbfl"
  outfl <- "testoutfl"
  on.exit({
    suppressWarnings(file.remove(file.path(ps[["wd"]], c(dbfl, outfl))))
  })
  run <- function(outfmt) {
    dbfl <- file.path(ps[["wd"]], "blast", dbfl)
    outfl <- file.path(ps[["wd"]], "blast", outfl)
    args <- c(
      "-query", dbfl, "-db", dbfl, "-outfmt", outfmt, "-dust",
      "no", "-strand", "plus", "-evalue", ps[["mxevl"]], "-out", outfl
    )
    cmdln(cmd = ps[["blstn"]], args = args, lgfl = dbfl, ps = ps)
  }
  sqs <- testsqs_gen(n = 3)
  blastdb_gen(sqs = sqs, dbfl = dbfl, ps = ps)
  outfmt_noquotes <- "6 qseqid sseqid pident length evalue qcovs qcovhsp"
  if (run(outfmt_noquotes) == 0) {
    return(outfmt_noquotes)
  }
  outfmt_quotes <- "\"6 qseqid sseqid pident length evalue qcovs qcovhsp\""
  if (run(outfmt_quotes) == 0) {
    return(outfmt_quotes)
  }
  error(ps = ps, "BLAST not working.")
}

# FILTER NOTES
# First filter out HSPs with lower than min.coverage
# Get query-subject pairs with too low coverage
# Remove both, query-subject and subject-query pair of low coverage
# hits!! Otherwise, we will end up with clusters with uneven sequence
# lengths!
# TODO does it make a difference to filter for qcovhsp ??
#' @name blast_filter
#' @title Filter BLAST results
#' @description Given a BLAST output, filters query-subject pairs
#' such that only HSPs with a coverage greater than \code{mncvrg}
#' (specified in the pipeline parameters) remain. Filters both:
#' query-subject and subject-query pairs, if one of the coverages is
#' insufficient. HSP coverage is obtained from the BLAST column
#' \code{qcovs}.
#' @param blast_res BLAST results
#' @template ps
#' @template lvl
#' @return data.frame blast res
#' @family run-private
blast_filter <- function(blast_res, ps, lvl = 3) {
  pull <- blast_res[["qcovs"]] < ps[["mncvrg"]]
  if (any(pull)) {
    # drop all with < mncvrg
    # ensure both qry-sbj and sbj-qry are dropped
    qsids <- blast_res[pull, c("query.id", "subject.id")]
    pull <- (blast_res[["query.id"]] %in% qsids[["subject.id"]] &
      blast_res[["subject.id"]] %in% qsids[["query.id"]]) | pull
  }
  ndrp <- sum(pull)
  info(
    lvl = lvl, ps = ps, "Removed [", ndrp, "/",
    nrow(blast_res), "] hits due to insufficient coverage"
  )
  if (sum(pull) > 0) {
    blast_res <- blast_res[!pull, ]
  }
  blast_res
}
