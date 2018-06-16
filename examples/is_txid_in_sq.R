data(tinamous)
sid <- tinamous@sids[[1]]
sq <- tinamous[[sid]]
txid <- sq@txid
# expect true
is_txid_in_sq(phylota = tinamous, txid = txid, sid = sid)
