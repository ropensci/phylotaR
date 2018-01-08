## ----setup, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  library(phylotaR)
#  wd <- '[YOUR PATH TO AOTUS FOLDER]'
#  ncbi_dr <- '[YOUR PATH TO NCBI BLAST TOOLS]'
#  txid <- 9504
#  setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr)

## ----running, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  run(wd=wd)

## ----restarting, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  reset(wd=wd, stage='cluster')
#  restart(wd=wd)

## ----selection, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  clstrs_obj <- genClstrsObj(wd)
#  clstr_id <- getBestClstrs(clstrs_obj, n=1)
#  writeSeqs(clstrs_obj, clstr_id, wd) # writes sequences to [clstr_id] + .fasta

## ----testing, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  library(ape)
#  system(paste0('mafft --auto ', clstr_id, '.fasta > alignment.fasta'))
#  system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n aotus -s alignment.fasta'))
#  tree <- read.tree(file='RAxML_bestTree.aotus')
#  plot(tree, no.margin=TRUE) # check if the same tax IDs are grouped

