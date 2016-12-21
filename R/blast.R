require('plyr')

make.blast.db <- function(seqs, dbfile='blastdb.fa', dir='.') {

    if (length(seqs) < 2) stop('Need more than 2 sequences')

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
    ## DUST filtering is enabled by default in blastn, disable
    ## Also only allow same-strand matches

    cmd <- paste('blastn -query', dbname, '-db', dbname, '-outfmt 6 -dust no -strand plus -evalue', evalue.cutoff, '-out', outfile)
    system(cmd) == 0 || stop('Command ', cmd, 'did not return 0')
    if (! file.exists(outfile)) {
        stop('Command ', cmd, 'did not produce output file ', outfile)
    }
    blast.results <- read.table(outfile)
    colnames(blast.results) <- c('query.id', 'subject.id', 'identity', 'alignment.length',
                                 'mismatches', 'gap.opens', 'q.start', 'q.end', 's.start',
                                 's.end', 'evalue', 'bit.score')
    return (blast.results)
}

filter.blast.results <- function(blast.results, seqs, min.coverage=0.51) {

    ## collapse HSPs such that we end up with unique query-subject pairs
    result.subset <- ddply(blast.results, c("query.id", "subject.id"), function(x)colSums(x['alignment.length']))

    ## adjust query and subject lengths for collapsed rows
    result.subset['query.length'] <- sapply(result.subset$query.id, function(id){seqs[[as.character(id)]]$length})
    result.subset['subject.length'] <- sapply(result.subset$subject.id, function(id){seqs[[as.character(id)]]$length})

    ## calculate how much overlap there is between hits
    coverages <- apply(result.subset, 1, function(x)x['alignment.length'] / max(x['query.length'], x['subject.length']))
    num.discarded.hits <- sum(coverages < min.coverage)
    cat("Discarding ", num.discarded.hits, " BLAST hits due to insufficient coverage\n")

    ## keep only the gis for which there is enough overlapping hits
    result.subset <- result.subset[which(coverages >= min.coverage),]

    return (result.subset)
}


require('RCurl')

url <- 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
program <- 'blastn'
database <- 'nt'
sequence <- URLencode('https://raw.githubusercontent.com/naturalis/supersmart/master/t/testdata/aln.fa')

rest.query <- 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastn&DATABASE=nt&QUERY=%3EFelis_margarita%0ATTGAGTACATGAATTGCACTTGGAACAGCAGCTCTGAGCCCCATCCCACCAACCTGACTCTGCACTATTG%0AGTATGAGAAGGGAAGAGGGGAGCACAGGGGAGTAAGAGGAGGTGGGCTGGATGGGACTTCGTGGGGACCA%0AAGAAAGAGGGTAGCCAGCATCCCAGCCTCCCTACCATTTTCTCATGGGGTAAGTCATACGTCAGTTCGGA%0AGATGAGGCTGGGCTGTCTTATCTGTAGTCCCCAAGTTTATACCACTGTTCCCCTTCCTCTCAACCCTTCT%0ACTAGGTACAAGAACTCCAATGATGATAGAGTCCAGGAG%0A'

result <- getURL(rest.query)

rid <- sub(".*RID = (.*?)\n.*", "\\1", result)
rtoe <- sub(".*RTOE = (.*?)\n.*", "\\1", result)

cat("Estimated waiting time: ", rtoe, " s", "\n")
while(TRUE) {
    Sys.sleep(1)

    query <- paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?ORGANISM=9685&CMD=Get&FORMAT_OBJECT=SearchInfo&RID=", rid)
    r <- getURL(query)

    wait <- grepl("WAITING", r)
    fail <- grepl("FAILED", r)
    unkn <- grepl("UNKNOWN", r)
    ready <- grepl("READY", r)

    cat("Waiting : ", wait, " fail : ", fail, " unknown : ", unkn, " ready : ", ready, "\n")

    if (ready == 1) {
        break
    }
}

getquery <-  paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=text&RID=", rid);
r <- getURL(getquery)

txt <- sub(".*?Value\n\n(.*?)\n\n.*", "\\1", r)
tab <- read.table(text=txt)

