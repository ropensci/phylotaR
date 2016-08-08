## script to extend the table we get with getnodes() of nodes of the
## NCBI taxonomy table. We add a column seqs, giving the number of genbank
## sequence records for each taxon

source('ncbi-tax.R')

## adds number of sequences for each taxon to table
add.numseqs <- function(taxid, nodes) {
    total <- nrow(nodes)
    num.seqs <- 0    
    for (child in children(taxid, taxdir)) {        
        ## 1 is its own parents, therefore omit
        if (child==1) {
            next
        }
        if (getrank(child, taxdir, nodes) %in% c('species', 'subspecies','varietas', 'forma')){
            child.num.seqs <- num.seqs.for.taxid(child)
            num.seqs <- num.seqs + child.num.seqs
            nodes[which(nodes$id==child),'seqs'] <<- child.num.seqs
            counter <<- counter + 1
            cat("Wrote nof seqs for (terminal) node # ", counter, "\n")
        }
        else {
            num.seqs <- num.seqs + add.numseqs(child, nodes)
        }
    }
    nodes[which(nodes$id==taxid),'seqs'] <<- num.seqs
    counter <<- counter + 1
    cat("Wrote nof seqs for (internal) node # ", counter, "\n")
    return(num.seqs)
}

counter <- 0
nodes <- ncbi.nodes
nodes$seqs <- 0

## add number of sequences for each node for root taxon
metazoa <- 33208
viridiplantae <- 33090
root.taxon <- metazoa

num.seqs <- add.numseqs(root.taxon, nodes)
save(nodes, file=paste0("nodes-seqnums-", root.taxon, '.Rdata'))





