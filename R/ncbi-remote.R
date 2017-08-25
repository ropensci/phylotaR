## functions for remote access to NCBI databases: Genbank and NCBI taxonomy
require('rentrez')

#' Provides a safer way to use 'rentrez' queries by retrying
#' a query in case of no success.
#' @param func A query function from package 'rentrez'
#' @param args A list with named arguments for 'func'
#' @param max.retry Number of tries to query remote server before givin gup
#' @return value returned by func or NULL if unsuccesful
#' @examples
#' search <- .safe.ncbi.query(entrez_search, list(db='nuccore', term="txid 1341354[Organism:exp]"))
.safe.ncbi.query <- function(func, args, max.retry=1000) {
    result <- NULL
    for (i in 1:max.retry) {
        ## Catch possible error
        query <- try(do.call(func, args))
        if( !is(query, "try-error") ) {
            result <- query
            break
        }
        else {
            ## retry if an error was thrown
            cat("Retry #", i, "calling rentrez function \n")
            print(args)
        }
    }
    return (result)
}

#' Given a taxon ID, get a vector with the Ganbank GI entries for this taxid
#' from the remote NCBI server.
#' @param taxid NCBI taxon identifier
#' @param max.len Maximum length of sequence to be included in GI list
#' @return vector with GIs
#' @examples
#' gis <- .gis.for.taxid(10149) ## get gis for capybara
.gis.for.taxid <- function(taxid, max.len=25000) {
    args <- list(db='nucleotide',
                 term=paste0('txid', taxid, '[Organism:noexp]', '1:', max.len, '[SLEN]'),
                 use_history=T,
                 retmax=1e9-1)
    search <- .safe.ncbi.query(entrez_search, args)

    return(search$ids)
}

#' Given a taxon ID, quieries NCBI for sequences that match the taxid.
#' @param taxid NCBI taxon identifier
#' @param max.len Maximum length of sequence to be included in GI list
#' @return list of lists containing sequence objects
#' @examples
#' seqs <- .seqs.for.taxid(10149) ## get gis for capybara
.seqs.for.taxid <- function(taxid, max.len=25000) {
    args <- list(db='nucleotide',
                 term=paste0('txid', taxid, '[Organism:noexp]', '1:', max.len, '[SLEN]'),
                 use_history=T,
                 retmax=1e9-1)
    search <- .safe.ncbi.query(entrez_search, args)
    gis <- search$ids;
    if (length(gis) < 1) {
        return(list())
    }
    cat("Going to retrieve ", length(gis), " sequences for taxid ", taxid, "\n")
    allseqs <- numeric()
    ## Fetch sequences in increments
    increment <- 500
    for (seq.idx in seq(0, length(gis), increment)) {
        ## Get FASTA strings for IDs in the specified segment
        cat("Retreiving seqs ", seq.idx+1, " to ",
            ifelse(seq.idx+increment<length(gis), seq.idx+increment, length(gis)),
            " for taxid ", taxid, "\n");
        seqargs <- list(db="nuccore",
                        rettype="fasta",
                        web_history=search$web_history,
                        retmax=increment,
                        retstart=seq.idx)
        res <- .safe.ncbi.query(entrez_fetch, seqargs)
        ## split fasta string on '>'
        seqstrs <- unlist(strsplit(res, "(?<=[^>])(?=\n>)", perl=T))
        ## clear defline
        seqstrs <- gsub("^\n?>.*?\\n", "", seqstrs)
        ## clear newlines
        seqstrs <- gsub("\n", "", seqstrs, fixed=TRUE)

        ## Get summary object to retreive taxid, accession, and other parameters for seqs
        summaries <- .safe.ncbi.query(entrez_summary, list(db='nucleotide',  web_history=search$web_history,retmax=increment, retstart=seq.idx))
        ## turn to list when only one result
        if (length(seqstrs)==1) {
            summaries <- list(summaries)
        }

        seqs <- lapply(1:length(seqstrs), function(i) {
            sm <- summaries[[i]]
            se <- seqstrs[[i]]
            seq <- list(gi=sm$uid,
                        ti=sm$taxid,
                        acc=sm$caption,
                        acc_vers=sm$accessionversion,
                        length=nchar(se),
                        ##TODO: Division not in summary object, how to get?
                        division=NA,
                        acc_date=sm$createdate,
                        ##TODO: GBrel not in summary object, how to get?
                        gbrel=NA,
                        def=sm$title,
                        seq=se)
            seq
        })
        cat("Done retreiving", length(seqs), "(", seq.idx+1, "to",
            ifelse(seq.idx+increment<length(gis), seq.idx+increment, length(gis)),
            ") seqs for taxid", taxid, "\n")

        allseqs <- c(allseqs, seqs)
    }
    cat("Finished retreiving", length(allseqs), "sequences for taxid", taxid, "\n")
    names(allseqs) <- gis
    return(allseqs)
}

#' Given a NCBI taxon ID, retreive the number of sequences for that taxon
#' @param taxid NCBI taxon identifier
#' @param max.len Maximum length of sequence to be included in count
#' @return numeric
.num.seqs.for.taxid <- function(taxid, max.len=25000) {
    args <- list(db='nucleotide',
                 term=paste0('txid', taxid, '[Organism:noexp]', '1:', max.len, '[SLEN]'),
                 use_history=T,
                 retmax=1e9-1)
    search <- .safe.ncbi.query(entrez_search, args)
    return(search$count)
}
