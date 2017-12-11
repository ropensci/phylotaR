library(phylotaR)
wd <- '/home/dom/Desktop/demo'
ncbi_execs <- setUpNcbiTools(d='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
txid <- 0000
setUpPrmtrs(wd=wd, txid=txid,
            ncbi_execs=ncbi_execs)
# Generate taxonomic 'nodes'
runTaxise(wd)
# Download sequences
runDownload(wd)
# Generate clusters
runClusters(wd)


