data('aotus')
txdct <- aotus@txdct
# this is a TaxDict object
# it contains taxonomic information, including records and tree
show(txdct)
# you can access its different data slots with @
txdct@txids  # taxonomic IDs
txdct@recs   # taxonomic records environment
txdct@txtr   # taxonomic tree
txdct@prnt   # MRCA
# access any record through the records environment
txdct@recs[[txdct@txids[[1]]]]
# for interacting with the taxonomic tree, see the treeman package
summary(txdct@txtr)
