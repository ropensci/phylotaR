library(phylotaR)
tdpth <- ''   # Provide the path to the NCBI folder, if excluded program will download automatically
wd <- ''      # Provide file path to working dir where output should be stored
ncbi_dr <- '' # Provide filepath to NCBI blast tools, if tools are on syspath, leave as is.
txid <- 9504  # NCBI txid of taxon choice
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, tdpth=tdpth)
run(wd=wd, nstages=3)
# reset(wd=wd, stage='cluster')  # pipeline can be reset and restarted
# restart(wd=wd, nstages=3)

# Early development of clusters object tools
clstrs_obj <- genClstrsObj(wd)
clstr_id <- getBestClstrs(clstrs_obj, n=1)
writeSeqs(clstrs_obj, clstr_id, wd) # writes sequences to [clstr_id] + .fasta

# if you want, run quick mafft and raxml
system(paste0('mafft --auto ', clstr_id, '.fasta > alignment.fasta'))
system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n aotus -s alignment.fasta'))
library(ape)
tree <- read.tree(file='RAxML_bestTree.aotus')
plot(tree, no.margin=TRUE) # check if the same taxa are grouped
