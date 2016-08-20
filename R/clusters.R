
cluster.seqs <- function(blast.results, informative=T) {
    g <- graph.data.frame(filtered.blast.results[,c("query.id", "subject.id")], directed=F)
    clusters <- clusters(g)

    ## filter for phylogenetically informative clusters
    if (informative) {
        clusters <- clusters$membership[which(clusters$membership %in% which (clusters$csize > 2))]
    }
    ## we will return a list, one entry with sequence IDs for each cluster
    cluster.list <- lapply(unique(informative.clusters), function(x)sort(names(informative.clusters)[which(informative.clusters==x)]))
    names(cluster.list) <- paste('cluster', 1:length(cluster.list), sep='')
    return (cluster.list)
}
