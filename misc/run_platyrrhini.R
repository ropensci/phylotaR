# 24 Nov 2017
# Current status: .dmp to manageable IDs
# I'm testing the pipeline with monkeys
# I like monkeys

# Set up
# Before running the pipeline you need to
#  - install NCBI blast tools
#  - look the path to the blast bin
#  - download NCBI taxonomy nodes.dmp and names.dmp
#  - place the files in a folder called 'NCBI' in the wd
# I'm working on R ways to automate these steps

# Set working directory
wd <- '/home/dom/Desktop/testing_phylotaR'
# Test if NCBI tools are present
ncbi_execs <- setUpNcbiTools(d='/home/dom/Programs/ncbi-blast-2.7.1+/bin')
# Select a txid, here Platyrrhini
txid <- 9479
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


