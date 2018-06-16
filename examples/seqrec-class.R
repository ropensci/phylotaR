data('aotus')
seqrec <- aotus@sqs@sqs[[1]]
# this is a SeqRec object
# it contains sequence records
show(seqrec)
# you can access its different data slots with @
seqrec@id       # sequence ID, accession + feature location
seqrec@nm       # feature name, '' if none
seqrec@accssn   # accession
seqrec@vrsn     # accession version
seqrec@url      # NCBI GenBank URL
seqrec@txid     # Taxonomic ID
seqrec@orgnsm   # free-text organism name
seqrec@sq       # sequence, in raw format
seqrec@dfln     # sequence definition
seqrec@ml_typ   # molecule type
seqrec@rec_typ  # whole record or feature
seqrec@nncltds  # sequence length
seqrec@nambgs   # number of non-ATCGs
seqrec@pambgs   # proportion of non-ATCGs
seqrec@gcr      # GC-ratio
seqrec@age      # days since being added to GenBank
# get the sequence like so....
(rawToChar(seqrec@sq))
