data('aotus')
taxrec <- aotus@txdct@recs[[aotus@txdct@txids[[1]]]]
# this is a TaxRec object
# it contains NCBI's taxonomic information for a single node
show(taxrec)
# you can access its different data slots with @
taxrec@id    # taxonomic ID
taxrec@scnm  # scientific name
taxrec@cmnm  # common name, '' if none
taxrec@rnk   # rank
taxrec@lng   # lineage information: list of IDs and ranks
taxrec@prnt  # parent ID
