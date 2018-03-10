# LIBS
library(phylotaR)

# LOAD USER-DEFINED PATHS
# set these paths for your own system
wd <- file.path(getwd(), 'demos', 'platyrrhini')
ncbi_dr <- readLines(file.path('demos', 'ncbi_dr.txt'))

# CREATE FOLDER, DELETE IF ALREADY EXISTS
if(file.exists(wd)) {
  unlink(wd, recursive=TRUE)
}
dir.create(wd)

# RUN PIPELINE
txid <- 9479
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, v=TRUE)
run(wd=wd)