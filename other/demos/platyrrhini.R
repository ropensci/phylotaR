library(phylotaR)
tdpth <- '/home/dom/Desktop/taxdump.tar.gz'
wd <- '/home/dom/Desktop/testing_phylotaR/platy2'
if(file.exists(wd)) {
  unlink(wd, recursive=TRUE)
}
dir.create(wd)
ncbi_dr <- '/home/dom/Programs/ncbi-blast-2.7.1+/bin'
txid <- 9479
file.create(file.path(wd, 'log.txt'))
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, tdpth=tdpth)
run(wd=wd, nstages=3)

# IDENTIFY BEST CLUSTER


# WRITE OUT AS SEQUENCES


# RUN MAFFT AND RAXML
