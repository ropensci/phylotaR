data('aotus')
# this is a Phylota object
# it contains cluster, sequence and taxonomic information from a phylotaR run
show(aotus)
# you can access its different data slots with @
aotus@cids   # cluster IDs
aotus@sids   # sequence IDs
aotus@txids  # taxonomic IDs
aotus@clstrs # clusters archive
aotus@sqs    # sequence archive
aotus@txdct  # taxonomic dictionary
# see all of the available slots
(slotNames(aotus))
# access different sequences and clusters with [[
(aotus[['0']])              # cluster record 0
(aotus[[aotus@sids[[1]]]])  # first sequence record
# get a summary of the whole object
(summary(aotus))
# the above generates a data.frame with information on each cluster:
# ID - unique id in the object
# Type - cluster type
# Seed - most connected sequence
# Parent - MRCA of all represented taxa
# N_taxa - number of NCBI recognised taxa
# N_seqs - number of sequences
# Med_sql - median sequence length
# MAD - Maximum alignment density, values close to 1 indicate all sequences in
#  the cluster have a similar length.
# Definition - most common words (and frequency) in sequence definitions
# Feature - most common feature name (and frequency)
