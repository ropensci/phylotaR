require('data.table')

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

blast.all.vs.all <- function(dbname='blastdb.fa', evalue.cutoff=1.0e-10, outfile='blastout.txt', dir='.') {

    outfile <- file.path(dir, outfile)
    dbname <- file.path(dir, dbname)

    ## Determine outformat. TODO: We don't really need all these columns...
    outfmt <- "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'"
    
    ## DUST filtering is enabled by default in blastn, disable
    ## Also only allow same-strand matches
        
    cmd <- paste('blastn -query', dbname, '-db', dbname, '-outfmt', outfmt, '-dust no -strand plus -evalue', evalue.cutoff, '-out', outfile)
    system(cmd) == 0 || stop('Command ', cmd, 'did not return 0')
    if (! file.exists(outfile)) {
        stop('Command ', cmd, 'did not produce output file ', outfile)
    }
    blast.results <- read.table(outfile)
    colnames(blast.results) <- c('query.id', 'subject.id', 'identity', 'alignment.length',
                                 'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start',
                                 's.end', 'evalue', 'bit.score', 'qcovs', 'qcovhsp')
    return (blast.results)
}

filter.blast.results <- function(blast.results, seqs, min.coverage=0.51) {

    ## recover()

    ## ( result.subset$q.end - result.subset$q.start) / result.subset$query.length * 100
    ## ( result.subset$s.end - result.subset$s.start) / result.subset$subject.length * 100
    
    ## collapse HSPs such that we end up with unique query-subject pairs
    ##result.subset <- ddply(blast.results, c("query.id", "subject.id"), function(x)colSums(x['alignment.length']))
    ##result.subset <- setDT(blast.results)[, .(alignment.length = sum(alignment.length)), .(query.id, subject.id)]

    ## first filter out HSPs with lower than min.coverage
    blast.dt <- data.table(blast.results)
   ## blast.dt <- blast.dt[qcovhsp/100 > min.coverage,]
    
    ## query-subject pairs with min coverage:
    blast.dt[qcovs/100 < min.coverage]    
    
    ## remove both, query-subject and subject-query pair!!
    ## Otherwise, we will end up with clusters with uneven sequences!!

    
    result.subset <- setkey(blast.dt[qcovs/100 > min.coverage,], query.id, subject.id)[, .(alignment.length = sum(alignment.length)), .(query.id, subject.id)]
    result.subset2 <- setkey(blast.dt[qcovhsp/100 > min.coverage,], query.id, subject.id)[, .(alignment.length = sum(alignment.length)), .(query.id, subject.id)]
    
    ##result.subset <- setkey(setDT(blast.results), query.id, subject.id)[, .(alignment.length = sum(alignment.length)), .(query.id, subject.id)]
    result.subset <- as.data.frame(result.subset)
    result.subset2 <- as.data.frame(result.subset2)
    ## adjust query and subject lengths for collapsed rows
    ## result.subset['query.length'] <- sapply(result.subset$query.id, function(id){seqs[[as.character(id)]]$length})
    ## result.subset['subject.length'] <- sapply(result.subset$subject.id, function(id){seqs[[as.character(id)]]$length})
    seqs.dt <- rbindlist(seqs)
    result.subset['query.length'] <- seqs.dt[match(result.subset$query.id, seqs.dt$gi),length]
    result.subset['subject.length'] <- seqs.dt[match(result.subset$subject.id, seqs.dt$gi),length]

    result.subset2['query.length'] <- seqs.dt[match(result.subset2$query.id, seqs.dt$gi),length]
    result.subset2['subject.length'] <- seqs.dt[match(result.subset2$subject.id, seqs.dt$gi),length]

    g <- graph.data.frame(result.subset[,c("query.id", "subject.id")], directed=F)
    g2 <- graph.data.frame(result.subset2[,c("query.id", "subject.id")], directed=F)
    
    ##hsp
    
    
##    seqs[as.character(result.subset$query.id)]
    
    ##result.subset['query.length'] <- unname(sapply(seqs[as.character(result.subset$query.id)], function(x)x$length))
    ##result.subset['subject.length'] <- unname(sapply(seqs[as.character(result.subset$subject.id)], function(x)x$length))


    ## calculate how much overlap there is between hits
    ##coverages <- apply(result.subset, 1, function(x)x['alignment.length'] / max(x['query.length'], x['subject.length']))

    ## we will discard all hits where the overlapping region is smaller
    ## than 0.51 (default) times of the length of both query and hit

    coverages <- result.subset[,'alignment.length'] /  pmax(result.subset[,'query.length'], result.subset[,'subject.length'])
    num.discarded.hits <- sum(coverages < min.coverage)
    cat("Discarding ", num.discarded.hits, " BLAST hits due to insufficient coverage\n")

    ## keep only the gis for which there is enough overlapping hits
    result.subset <- result.subset[which(coverages >= min.coverage),]

    return (result.subset)
}


