## ----setup, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  library(phylotaR)
#  wd <- '[YOUR PATH TO AOTUS FOLDER]'
#  ncbi_dr <- '[YOUR PATH TO NCBI BLAST TOOLS]'
#  txid <- 9504
#  setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, v=TRUE)

## ----running, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  run(wd=wd)

## ----restarting, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  reset(wd=wd, stage='cluster')
#  restart(wd=wd)

## ----selection1, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  all_clustrs <- read_phylota(wd)
#  cids <- all_clustrs@cids
#  n_taxa <- get_ntaxa(phylota=all_clustrs, cid=cids)

## ----selection2, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  keep <- cids[n_taxa > 6]
#  selected <- drop_cls(phylota=all_clustrs, cid=keep)
#  smmry <- summary(selected)

## ----selection3, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  cid <- smmry[2, 'ID']
#  # get the cluster record
#  cluster_record <- selected@cls[[cid]]
#  # use the seq. IDs to get the sequence records
#  seq_records <- selected@sqs[cluster_record@sids]
#  # extract a single record
#  seq_record <- seq_records[[seq_records@ids[[1]]]]
#  # get the sequence
#  seq <- rawToChar(seq_record@sq)

## ----selection4, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  # choose best sequence per species
#  reduced <- drop_by_rank(phylota=selected, rnk='species', n=1)
#  # get txids at the species level for each sequence
#  txids <- get_txids(phylota=reduced, cid=cid, rnk='species')
#  # look up name for txids
#  scientific_names <- get_tx_slot(phylota=reduced, txid=txids, slt_nm='scnm')
#  # clean the names
#  scientific_names <- gsub('\\.', '', scientific_names)
#  scientific_names <- gsub('\\s+', '_', scientific_names)
#  # look up sequence IDs for our chosen cluster
#  sids <- reduced@cls[[cid]]@sids
#  # write out
#  write_sqs(phylota=reduced, sid=sids, outfile='cytb.fasta',
#            sq_nm=scientific_names)

## ----testing, eval=FALSE, message=FALSE, warning=FALSE, include=TRUE, paged.print=FALSE----
#  library(ape)
#  system('mafft --auto cytb.fasta > alignment.fasta')
#  system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n aotus -s alignment.fasta'))
#  tree <- read.tree(file='RAxML_bestTree.aotus')
#  plot(tree, no.margin=TRUE, type='unrooted')

