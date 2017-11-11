
## TODO: maybe need unique filename, maybe add if it is a direct or subtree cluster to calculate
make.blast.db <- function(seqs, dbfile='blastdb.fa', dir='.') {

    if (length(seqs) < 2) stop('Need more than 2 sequences')
    cat("Making blast database for ", length(seqs), " sequences\n")

    file <- file.path(dir, dbfile)
    file.create(file)

    ## Write user-supplied gis to file
    for (gi in names(seqs)) {
        write(paste0("> ", gi, "\n", seqs[[as.character(gi)]]$seq), file=file, append=T)
    }
    cmd <- paste('makeblastdb -in', file, '-dbtype nucl')
    system(cmd) == 0 || stop('Command ', cmd, 'did not return 0')

    ## Check if files produced by makeblastdb command are present
    extensions <- c('nhr', 'nin', 'nsq')
    fnames <- sapply(extensions, function(e)paste0(file, '.', e))
    if ( ! all (sapply(fnames, file.exists))) {
        stop('Command ', cmd, ' did not produce output files ', paste(fnames))
    }
    return (file)
}

## TODO: maybe need unique filename, maybe add if it is a direct or subtree cluster to calculate
blast.all.vs.all <- function(dbname='blastdb.fa', evalue.cutoff=1.0e-10, outfile='blastout.txt', dir='.') {

    outfile <- file.path(dir, outfile)
    dbname <- file.path(dir, dbname)

    ## Determine outformat. TODO: We don't really need all these columns...
    outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'"

    ## DUST filtering is enabled by default in blastn, disable
    ## Also only allow same-strand matches
    cmd <- paste('blastn -query', dbname, '-db', dbname, '-outfmt', outfmt, '-dust no -strand plus -evalue', evalue.cutoff, '-out', outfile)
    system(cmd) == 0 || stop('Command ', cmd, 'did not return 0')

    ## catch no-result case
    if (! file.exists(outfile)) {
        cat("No BAST output, returning NULL\n")
        return(NULL)
    }
    ## TODO: cannot include this above because it can produce NA
    if (file.info(outfile)$size==0) {
        cat("No BAST output, returning NULL\n")
        return(NULL)
    }

    blast.results <- read.table(outfile)
    colnames(blast.results) <- c('query.id', 'subject.id', 'identity', 'alignment.length',
                                 'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start',
                                 's.end', 'evalue', 'bit.score', 'qcovs', 'qcovhsp')
    return (blast.results)
}

filter.blast.results <- function(blast.results, seqs, min.coverage=0.51) {

    blast.dt <- data.table(blast.results)

    ## First filter out HSPs with lower than min.coverage
    ## Get query-subject pairs with too low coverage
    ## Remove both, query-subject and subject-query pair of low coverage hits!!
    ## Otherwise, we will end up with clusters with uneven sequence lengths!
    ## TODO does it make a difference to filter for qcovhsp ??
    remove.dt <- blast.dt[qcovs/100 < min.coverage, .(query.id, subject.id)][, rbind(.SD, rev(.SD), use.names = FALSE)]
    blast.dt.filtered <- blast.dt[!remove.dt, on=names(remove.dt)]
    cat("Removed", nrow(blast.dt) - nrow(blast.dt.filtered), "of the", nrow(blast.dt), "BLAST hits due to insufficient coverage\n")

    ## XXX Don't need to collapse when doing single-linkeage clustering!!
    ## collapse HSPs such that we end up with unique query-subject pairs
    ## TODO The xx is a hack.. how do I collapse by two column values and only get the column back??
    ## result.subset <- setkey(blast.dt.filtered, query.id, subject.id)[, .(xx=0),.(query.id, subject.id)]
    result <- as.data.frame(blast.dt.filtered[,.(query.id, subject.id)])

    return (result)
}


