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

reset(wd=wd, stage='cluster')
restart(wd=wd, nstages=3)

# dev of cluster obj
clstrs_obj <- genClstrsObj(wd)
clstr_id <- getBestClstrs(clstrs_obj, n=1)
writeSeqs(clstrs_obj, clstr_id, wd) # writes sequences to [clstr_id] + .fasta

# run quick mafft and raxml
system(paste0('mafft --auto ', clstr_id, '.fasta > alignment.fasta'))
system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n aotus -s alignment.fasta'))
library(ape)
tree <- read.tree(file='RAxML_bestTree.aotus')
plot(tree, no.margin=TRUE) # check if the same taxa are grouped
