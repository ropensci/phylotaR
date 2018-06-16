data('aotus')
clstrrec <- aotus@clstrs@clstrs[[1]]
# this is a ClstrRec object
# it contains cluster information
show(clstrrec)
# you can access its different data slots with @
clstrrec@id     # cluster id
clstrrec@sids   # sequence IDs
clstrrec@nsqs   # number of sequences
clstrrec@txids  # taxonomic IDs of sequences
clstrrec@ntx    # number unique taxonomic IDs
clstrrec@typ    # cluster type: merged, subtree, direct or paraphyly
clstrrec@prnt   # MRCA of all taxa
clstrrec@seed   # most inter-connected sequence
