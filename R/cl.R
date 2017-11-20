# CIRCULAR REF
#source('blast.R')
#source('ncbi-remote.R')

## Schema in phylota database:
## CREATE TABLE "clusters_194" (
##  "ti_root" int(10)  DEFAULT NULL,
##  "ci" int(10)  DEFAULT NULL,
##  "cl_type" text  DEFAULT NULL,
##  "n_gi" int(10)  DEFAULT NULL,
##  "n_ti" int(10)  DEFAULT NULL,
##  "PI" tinyint(1) DEFAULT NULL,
##  "MinLength" int(10)  DEFAULT NULL,
##  "MaxLength" int(10)  DEFAULT NULL,
##  "MaxAlignDens" float DEFAULT NULL,
##  "ci_anc" int(10)  DEFAULT NULL,
##  "seed_gi" bigint(20)  DEFAULT NULL,
##  "Q" float DEFAULT NULL,
##  "TC" float DEFAULT NULL,
##  "clustalw_tree" longtext,
##  "muscle_tree" longtext,
##  "strict_tree" longtext,
##  "clustalw_res" float DEFAULT NULL,
##  "muscle_res" float DEFAULT NULL,
##  "strict_res" float DEFAULT NULL,
##  "ortho" tinyint(4) DEFAULT NULL,
##  "n_gen" int(10)  DEFAULT NULL,
##  "n_child" int(10)  DEFAULT NULL,
##  "muscle_tree_mid" longtext,
##  "muscle_tree_mid_chron" longtext,
##  "muscle_tree_aln_length" int(10)  DEFAULT NULL,
##  "PD_data" float DEFAULT NULL,
##  "PD_time" float DEFAULT NULL,
##  "brlen_time_max" float DEFAULT NULL,
##  "mp_score" int(10)  DEFAULT NULL,
##  "uid" int(11) NOT NULL ,
##  "uid_anc" int(10)  DEFAULT NULL,
##  "root_cluster" int(10)  DEFAULT NULL,
##  "lca" int(10)  DEFAULT NULL,
##  PRIMARY KEY ("uid")
##);

#' @name clusters.ci_gi.seqs.create
#' @title Clusters CI GI Sequences
#' @description TODO
#' @details
#' @export
#' @examples
#' # TODO
clusters.ci_gi.seqs.create <- function(root.taxon, nodesfile,
                                       files=list(clusters='clusters.tsv',
                                                  ci_gi='ci_gi.tsv',
                                                  seqs='seqs.tsv'),
                                       informative=FALSE,
                                       prmtrs=prmtrs) {
    ## load nodes table
    nodes <- read.table(nodesfile, header=TRUE)

    ## load clusters that are already calculated, if any
    clusters.file <- files[['clusters']]
    ci_gi.file <- files[['ci_gi']]
    seqs.file <- files[['seqs']]

    ## get all nodes that will be processed:
    ## all nodes for which all children contain < MAX.BLAST.SEQS sequences
    taxa.to.process <- vector()
    queue <- root.taxon
    while(length(queue) > 0) {
        current.taxon <- head(queue, 1)
        queue <- tail(queue, length(queue)-1)

        ## get number of sequences to determine if it is manageable to calculate clusters
        ## TODO: We should get the counts from the nodes table!!!
        cat("Counting species for taxon", current.taxon, "\n")
        num.nonmodel.seqs <- nodes[match(current.taxon, nodes$ti),
                                   'n_gi_sub_nonmodel']
        num.model.seqs <- nodes[match(current.taxon, nodes$ti),
                                'n_sp_model'] * prmtrs[['mdl_thrshld']]
        num.seqs <- num.nonmodel.seqs + num.model.seqs

        cat("Number of sequences for taxon", current.taxon, ": ",
            num.seqs, "\n")

        ## if sequence count is smaller than MAX.BLAST.SEQS,
        # add it to the nodes to process, otherwise get children
        if (num.seqs <= prmtrs[['mx_blst_sqs']]) {
            cat("Will process taxon", current.taxon, "\n")
            taxa.to.process <- c(taxa.to.process, current.taxon)
        }
        else {
            cat("Too many seqs to blast for taxid", current.taxon,
                "... retrieving children\n")
            queue <- c(queue, .children(current.taxon, nodes))
        }
    }

    foreach (i=seq_along(taxa.to.process), .verbose=TRUE) %dopar% {
    ##for (i in seq_along(taxa.to.process)) {
        taxid <- taxa.to.process[i]
        cat("Processing taxid ", taxid, " # ", i, " / ",
            length(taxa.to.process), "\n")

        ## Check if we can cache sequences for this taxon
        seqs <- .read.seqs.from.fasta(taxid, prmtrs[['sq_cch_dr']],
                                      prmtrs[['mdl_thrshld']])

        ## Get the sequences if they could not be read from file
        if (length(seqs)==0) {
            seqs <- .rec.retrieve.seqs(taxid, nodes,
                                       max.len=prmtrs[['mx_sq_lngth']],
                                       max.seqs.per.spec=prmtrs[['mdl_thrshld']])
            .write.seqs.to.fasta(seqs, taxid, prmtrs[['sq_cch_dr']],
                                 prmtrs[['mdl_thrshld']])
        }

        cat("Making sequence dataframe\n")
        seqdf <- do.call(rbind, lapply(seqs, as.data.frame))
        cat("Done making sequence dataframe\n")
        clusters <- cluster(taxid, nodes, seqs, informative=informative)
        cldf <- .make.cluster.entries(clusters)
        cigidf <- .make.ci_gi.entries(clusters)

        ## Write all data to file
        cat("Taxid", taxid, ": writing", nrow(cldf), "clusters,",
            nrow(seqdf), "sequences,", nrow(cigidf), "ci_gi entries to file\n")
        write.table(cldf, file=clusters.file, append=file.exists(clusters.file),
                    quote=FALSE, sep="\t", row.names=FALSE,
                    col.names=!file.exists(clusters.file))
        write.table(seqdf, file=seqs.file, append=file.exists(seqs.file),
                    quote=FALSE, sep="\t", row.names=FALSE,
                    col.names=!file.exists(seqs.file))
        write.table(cigidf, file=ci_gi.file, append=file.exists(ci_gi.file),
                    quote=FALSE, sep="\t", row.names=FALSE,
                    col.names=!file.exists(ci_gi.file))
        cat("Finished processing taxid ", taxid, " # ", i, " / ",
            length(taxa.to.process), "\n")
    }
}

## input : List of sequences
## filename is {taxon}-max-{threshold}.fa
.write.seqs.to.fasta <- function(seqs, taxon, dir, threshold=NA) {
    fastastr <- ""

    ## assemble fasta string
    for (s in seqs) {
        ## defline will consist of all attributes but 'seq' separated by "|".
        # Example:
        ## gi=1104556983|ti=61454|acc=KX265095|acc_vers=KX265095.1 ...
        attributes <- names(s)[!grepl("seq", names(s))]
        defline <- paste(">", paste(paste(attributes, s[attributes], sep="="),
                                    collapse="|"))
        ## append to fasta string
        fastastr <- paste0(fastastr, defline, "\n", s$seq, "\n")
    }
    ## write to file
    filename <- paste0(dir, "/", taxon, "-max-", threshold, ".fa")
    cat("Writing", length(seqs), "sequences for taxon", taxon, "to file",
        filename, "\n")
    if (! file.exists(dir)) {
        dir.create(dir)
    }
    cat(fastastr, file=filename)
}

.read.seqs.from.fasta <- function(taxon, dir, threshold=NA) {
    all.seqs <- list()

    ## get filename and read in fasta
    filename <- paste0(dir, "/", taxon, "-max-", threshold, ".fa")

    ## check if file exists
    if (! file.exists(filename)) {
        return(all.seqs)
    }

    cat("Reading sequences for taxid", taxon, "from file", filename, "\n")

    ## each line in file is an item in the vector
    vec <- scan(filename, what=character(), sep="\n")
    deflines <- vec[seq(1, length(vec), by=2)]
    seqstrs <- vec[seq(2, length(vec), by=2)]

    ## parse sequence info from deflines
    deflines <- sub("^> ", "", deflines)
    ll <- lapply(deflines, function(d)strsplit(d, "\\|")[[1]])

    attributes <- c("gi", "ti", "acc", "acc_vers", "length", "division",
                    "acc_date", "gbrel", "def")

    for (i in seq_along(ll)) {
        ## TODO: This is a bit dirty, works only if there is no "|"
        # character in defline!
        values <- unname(sapply(attributes,
                                function(x){ gsub(paste0(".*?", x, "=(.*?)\\|.*"),
                                                  "\\1", paste0(deflines[i], "|"))}))
        seq <- list()
        seq[attributes] <- values
        seq$seq <- seqstrs[[i]]
        all.seqs[[i]] <- seq
    }
    names(all.seqs) = lapply(all.seqs, `[[`, "gi")
    cat("Read", length(all.seqs), "sequences for taxid" , taxon,
        "from file", filename, "\n")
    return (all.seqs)
}

#' @name Cluster
#' @title Cluster
#' @description TODO
#' @details
#' @export
#' @examples
#' # TODO
cluster <- function(taxon, nodes, seqs, blast.results=NULL,
                    direct=FALSE, informative=FALSE) {
    ## list with clusters to be returned
    all.clusters <- list()

    rank <- as.character(nodes$rank[match(taxon, nodes$ti)])
    cat ("Processing taxid", taxon, "of rank", rank, "attempting to make",
         ifelse(direct, "direct", "subtree"), "clusters\n")

    ## Taxa can have "direct" sequence links up to some certain level, e.g. genus.
    # This means that no sequnces of the children
    ## of the taxon are included. These clusters will be stored separately
    ## as "node" clusters. In genbank, one searches for direct links using keyword
    # 'noexp' and 'exp' for subtree links.
    ## Below, we also calulate the direct ('node') clusters for the taxon of
    # interest
    if (! direct) {

        ## In NCBI, subtree search on terminals (taxa without children)
        # gives the direct sequences.
        ## Therefore, we have to exclude terminals from direct cluster calculation,
        # because their
        ## clusters are already calculated in the 'subtree' (direct=false) mode!
        if (length(.children(taxon, nodes)) > 0) {
            all.clusters <- c(all.clusters, cluster(taxon, nodes, seqs,
                                                    blast.results, TRUE,
                                                    informative))
        }
    }
    ## get GIs for taxon
    taxids <- c()
    if (direct) {
        taxids <- taxon
    }
    else {
        taxids <- .get.subtree.taxids(taxon, nodes)
    }
    all.tis <- sapply(seqs, '[[', 'ti')
    gis <- names(seqs[which(all.tis %in% taxids)])
    cat("Found", length(gis), ifelse(direct, "direct", "subtree"),
        "gis for taxid", taxon, "\n")

    if (length(gis) == 0) {
        cat("No", ifelse(direct, "direct", "subtree"), "sequences for taxid",
            taxon, ",cannot make clusters\n")
        return(all.clusters)
    }
    current.seqs <- seqs[gis[gis %in% names(seqs)]]

    ## Get blast results if not already calculated
    current.blast.results <- list()
    if (is.null(blast.results)) {
        cat("Current BLAST result is NULL from parent, performing BLAST\n")
        current.blast.results <- .get.blast.results(taxon, current.seqs)

        ## exit function if still no blast results
        ## TODO: Can we do this more elegantly?
        if (is.null(current.blast.results)) {
            cat("Current BLAST result is NULL after blasting\n")
            return(all.clusters)
        }
    }
    else {
        ##  reduce blast results such that include only the gis for
        # the current taxon
        cat("Using BLAST results from parent cluster\n")
        current.blast.results <- blast.results[which(
          blast.results$subject.id %in% gis),]
        current.blast.results <- current.blast.results[which(
          current.blast.results$query.id %in% gis),]
    }

    ## Sequence clusters are stored in a list of lists, named by gi
    ## Get top-level clusters
    cat("Clustering BLAST results for taxid", taxon, "\n")
    raw.clusters <- cluster.blast.results(current.blast.results,
                                          informative=informative)
    cat("Finished clustering BLAST results for taxid", taxon, "\n")

    ## make dataframe with fields as in PhyLoTa database
    clusters <- .add.cluster.info(raw.clusters, nodes, taxon, seqs,
                                  direct=direct)
    cat("Generated ", length(clusters), "clusters\n")
    all.clusters <- c(all.clusters, clusters)

    ## If we do not calculate a direct cluster here, iterate over
    # taxon's children to retrieve clusters
    if (! direct) {
        for (ch in .children(taxon, nodes)) {
            cat("Processing child taxon of ", taxon, " ", ch, "\n")

            ## calculate clusters for child taxon
            child.clusters <- cluster(ch, nodes, seqs, current.blast.results)
            ## two numbers can only be calculated if we have cluster
            # info on multiple taxonomic levels:
            ##  n_child, the number of child clusters, and ci_anc,
            # the parent cluster.
            ##  Since we are dealing with single-likeage clusters,
            # we can identify the parent cluster of a
            ##  cluster if at least one gi of both clusters is the same.
            child.clusters <- lapply(child.clusters, function(cc) {
                ## Get the indices of the parent clusters that contain
                #a gi of the child clusters
                ## If a gi is in two clusters (e.g. parent and grandparent),
                # take the lower one, by taking
                ## the higher index
                idx <- max(which(sapply(clusters, function(c) any(
                  cc$gis %in% c$gis))))
                ## set parent cluster to child cluster
                cc$ci_anc <- clusters[[idx]]$ci
                ## update child count of parent clusters
                clusters[[idx]]$n_child <- clusters[[idx]]$n_child + 1
                cc
            })
            ## append child clusters to result cluster list
            all.clusters <- c(all.clusters, child.clusters)
        }
    }
    return(all.clusters)
}



.rec.retrieve.seqs <- function(taxid, nodes, max.len=25000,
                               max.seqs.per.spec=10000) {
    cat("Attempting to retrieve sequences for taxid", taxid, "\n")
    seqs <- list()
    ##node.ranks <- c('species', 'subspecies', 'varietas', 'forma')

    ## get subtree counts. If that is smaller than max.seqs.per.spec,
    # we are done and these are the sequences
    subtree.count <- .num.seqs.for.taxid(taxid, direct=FALSE,
                                         max.len=max.len)
    if (subtree.count <= max.seqs.per.spec) {
        cat(subtree.count, "seqs for taxon", taxid, ", less than maximum of ",
            max.seqs.per.spec, " ")
        cat("sequences. Retreiving sequences for whole subtree\n")
        seqs <- .seqs.for.taxid(taxid, direct=FALSE, max.len=max.len,
                                max.seqs=max.seqs.per.spec)
        return (seqs)
    }

    ## First retrieve direct sequence of focal taxon,
     # then the ones of the children
    seqs <- c(seqs, .seqs.for.taxid(taxid, direct=TRUE,
                                    max.len=max.len, max.seqs=max.seqs.per.spec))
    for (ch in .children(taxid, nodes)) {
        seqs <- c(seqs, .rec.retrieve.seqs(ch, nodes, max.len,
                                           max.seqs.per.spec))
    }
    return(seqs)
}

## TODO: Can this function be deleted?
.rec.num.seqs <- function(taxid, nodes, max.len=25000,
                          max.seqs.per.spec=100000) {
    count <- 0

    ## get subtree counts. If that is smaller than max.seqs.per.spec,
    # we don't need to
    ##  traverse the taxonomy tree
    subtree.count <- .num.seqs.for.taxid(taxid, direct=FALSE,
                                         max.len=max.len)
    if (subtree.count <= max.seqs.per.spec) {
        return (subtree.count)
    }

    ## retrieve sequence counts
    c <- .num.seqs.for.taxid(taxid, direct=TRUE,
                             max.len=max.len)

    ## if exceeded maximum amount per spec (or lower),
    # return the count that will be retrieved
    if (c > max.seqs.per.spec) {
        c <- max.seqs.per.spec
    }
    count <- count + c

    for (ch in .children(taxid, nodes)) {
        count <- count + .rec.num.seqs(ch, nodes, max.len,
                                       max.seqs.per.spec)
    }
    return(count)
}

.get.blast.results <- function(taxon, seqs) {
    cat("BLAST all vs all for", length(seqs), "sequences\n")
    ## make unique name for BLAST files
    ## TODO: Change this to tempdir later
    dir = 'blast'
    if (! file.exists(dir)) {
        dir.create(dir)
    }
    dbfile <- paste0('taxon-', taxon, '-db.fa')
    outfile <- paste0('taxon-', taxon, '-blastout.txt')
    make.blast.db(seqs, dbfile=dbfile, dir=dir)
    blast.results <- blast.all.vs.all(dbfile,
                                      outfile=outfile, dir=dir)
    ## TODO: Not so elegant
    if (is.null(blast.results)) {
        return(NULL)
    }
    cat("Number of BLAST results ", nrow(blast.results), "\n")
    cat("Filtering BLAST results\n")
    filtered.blast.results <- filter.blast.results(blast.results, seqs)

    return(filtered.blast.results)
}

#' @name cluster.blast.results
#' @title Cluster BLAST Results
#' @description TODO
#' @details
#' @export
#' @examples
#' # TODO
cluster.blast.results <- function(blast.results, informative=TRUE) {
    g <- graph.data.frame(blast.results[,c("query.id", "subject.id")],
                          directed=FALSE)
    clusters <- clusters(g)

    ## filter for phylogenetically informative clusters
    if (informative) {
        clusters <- clusters$membership[which(
          clusters$membership %in% which (clusters$csize > 2))]
    } else {
        clusters <- clusters$membership
    }
    ## we will return a list, one entry with sequence IDs for each cluster
    cluster.list <- lapply(unique(clusters), function(x) {
        list(gis=sort(names(clusters)[which(clusters==x)]))
    })

    ## Get the seed gi, we will chose it to be the sequence
    # in the cluster that has
    ## the most hits with the other members in the cluster;
    # i.e. the most connected
    ## node in the graph
    degrees <- degree(g)
    ## get seed gis and as field to clusters
    cluster.list <- lapply(cluster.list, function(cl){
        idx <- order(degrees[cl$gis], decreasing=T)[1]
        ## index of most connected component
        cl$seed_gi <- cl$gis[idx]
        cl
    })
    return (cluster.list)
}

## input: List of clusters
.add.cluster.info <- function(clusters, nodes, taxid,
                              seqs, direct=FALSE) {
     ## we will take the column names in the database as names in the list
    for (i in 1:length(clusters)) {
        cl <- clusters[[i]]
        cl.seqs <- lapply(cl$gis, function(gi)seqs[[as.character(gi)]])
        cl$ti_root <- taxid
        cl$ci <- i-1
        cl$cl_type <- ifelse(direct, 'node', 'subtree')
        cl$n_gi <- length(cl$gis)
        ## all taxon ids for gis in cluster
        cl$tis <- sapply(cl.seqs, '[[', 'ti')
        cl$n_ti <- length(unique(cl$tis))
        ## sequence lengths
        l <- sapply(cl.seqs, '[[', 'length')
        cl$MinLength <- min(l)
        cl$MaxLength <- max(l)
        ## get number of genera
        cl$n_gen <- length(unique(sapply(cl$tis, .genus.for.taxid,
                                         nodes)))
        ## n_child and ci_anc will be set (or incremented) later,
        # when multiple hierarchies are calculated
        cl$n_child <- 0
        cl$ci_anc <- NA
        ## make a unique cluster id consisting of seed gi, taxon id,
        # cluster id, cluster type
        unique.id <- paste0(cl$seed_gi, '-', taxid, '-',
                            cl$ci, '-', cl$cl_type)
        cl$unique_id <- unique.id
        clusters[[i]] <- cl
    }
    return(clusters)
}

## returns a data frame with columns matchine the fields in
# PhyLoTA's cluster table
.make.cluster.entries <- function(clusters) {
    ## remove gi and id fields which won't be part of cluster table
    l <- lapply(clusters, function(cl)cl[which(!names(cl)%in%c(
      'gis', 'tis', 'unique_id'))])
    ## create data frame
    df <- do.call(rbind, lapply(l, as.data.frame))
    return(df)
}

.get.subtree.taxids <- function(taxid, nodes) {
    children <- .children(taxid, nodes)
    result <- c(taxid)
    for (ch in children) {
        result <- c(result, .get.subtree.taxids(ch, nodes))
    }
    return(result)
}

.children <- function(taxid, nodes) {
    nodes[nodes[,'ti_anc']==taxid,'ti']
}

.genus.for.taxid <- function(taxid, nodes) {
    nodes[nodes[,'ti_anc']==taxid,'ti_genus']
}

