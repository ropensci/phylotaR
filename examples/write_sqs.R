data('aotus')
# get sequences for a cluster and write out
random_cid <- sample(aotus@cids, 1)
sids <- aotus[[random_cid]]@sids
write_sqs(phylota = aotus, outfile = 'test.fasta', sq_nm = 'my_gene',
          sid = sids)
