data("tinamous")
# Plot clusters, size by n. sq, fill by n. tx
p <- plot_phylota_treemap(phylota = tinamous, cids = tinamous@cids,
                          area = 'nsq', fill = 'ntx')
print(p)
# Plot taxa, size by n. sq, fill by ncl
txids <- get_txids(tinamous, txids = tinamous@txids, rnk = 'genus')
txids <- txids[txids !=  '']
txids <- unique(txids)
txnms <- get_tx_slot(tinamous, txids, slt_nm = 'scnm')
p <- plot_phylota_treemap(phylota = tinamous, txids = txids, txnms = txnms,
                          area = 'nsq', fill = 'ncl')
print(p)
