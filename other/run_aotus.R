library(phylotaR)
tdpth <- '/home/dom/Desktop/taxdump.tar.gz'
wd <- '/home/dom/Desktop/testing_phylotaR/aotus'
if(file.exists(wd)) {
  unlink(wd, recursive=TRUE)
}
dir.create(wd)
ncbi_dr <- '/home/dom/Programs/ncbi-blast-2.7.1+/bin'
txid <- 9504
file.create(file.path(wd, 'log.txt'))
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, tdpth=tdpth)
run(wd=wd, nstages=3)
