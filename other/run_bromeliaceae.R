library(phylotaR)

# Set working directory
wd <- '/home/dom/Desktop/testing_phylotaR/brom'
# Test if NCBI tools are present
ncbi_execs <- setUpNcbiTools(d='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
# Select a txid, here Platyrrhini
txid <- 4613
setUpPrmtrs(wd=wd, txid=txid,
            ncbi_execs=ncbi_execs,
            mx_dscndnts=10000,
            tmout=100,
            mdl_thrshld=3000,
            mx_blst_sqs=10000,
            mx_sq_lngth=25000,
            verbose=TRUE,
            cores=2)
# Generate taxonomic 'nodes'
genTxNds(wd)
