# Set working directory
wd <- '/home/dom/Desktop/testing_phylotaR'
# FYI: requires NCBI taxonomy dmp to be in a folder
# called NCBI/
# Set up execs and parameters
ncbi_execs <- setUpNcbiTools(d='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
txid <- 9479
setUpPrmtrs(wd=wd, txid=txid,
            ncbi_execs=ncbi_execs,
            mdl_thrshld=3000,
            mx_blst_sqs=10000,
            mx_sq_lngth=25000,
            verbose=TRUE,
            cores=2)


