data('aotus')
seqarc <- aotus@sqs
# this is a SeqArc object
# it contains sequence records
show(seqarc)
# you can access its different data slots with @
seqarc@ids     # sequence IDs defined as accession + feature position
seqarc@nncltds # number of nucleotides of all sequences  
seqarc@nambgs  # number of ambiguous nucleotides of all sequences
seqarc@txids   # all the taxonomic IDs for all sequences
seqarc@sqs     # list of all SeqRecs
# access sequence records [[
(seqarc[[seqarc@ids[[1]]]])  # first sequence record
# generate new sequence archives with [
(seqarc[seqarc@ids[1:10]])  # first 10 sequences
