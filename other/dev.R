
devtools::load_all('~/Coding/phylotaR')

ps <- ldPrmtrs('~/Desktop/big_tree_of_birds/phylotar')



wd <- '/home/dom/Coding/phylotaR_demo/primates'
ps <- ldPrmtrs(wd)
ps[['wd']] <- wd
ps[['v']] <- TRUE
ps$mkblstdb <- aotus_ps$mkblstdb
ps$blstn <- aotus_ps$blstn
clstrClstrs(ps)

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")


