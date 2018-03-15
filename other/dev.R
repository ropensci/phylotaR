
devtools::load_all('~/Coding/phylotaR')

data(aotus)
phylota <- aotus
phylota <- drop_cls(phylota, phylota@cids[1:10])

cids <- phylota@cids
txids <- phylota@txids
txnms <- get_tx_slot(phylota, txids, slt_nm='scnm')
txnms <- sample(txnms)
cnms <- cids


plot_phylota_pa <- function(phylota, cids, txids,
                            cnms=cids, txnms=txids) {
  mkdata <- function(cid) {
    cl <- phylota@cls@cls[[which(cid == phylota@cids)]]
    value <- factor(as.numeric(txids %in% cl@txids))
    data.frame(txid=as.character(txids),
               cid=cid, value=value)
  }
  # gen p_data
  p_data <- lapply(cids, mkdata)
  p_data <- do.call('rbind', p_data)
  p_data[['cnm']] <- 
    cnms[match(p_data[['cid']], cids)]
  p_data[['txnm']] <- 
    txnms[match(p_data[['txid']], txids)]
  # reorder
  p_data[['cnm']] <- factor(p_data[['cnm']],
                            levels=cnms,
                            ordered=TRUE)
  p_data[['txnm']] <- factor(p_data[['txnm']],
                             levels=txnms,
                             ordered=TRUE)
  # plot
  ggplot2::ggplot(p_data, ggplot2::aes(cnm, txnm)) +
    ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::scale_fill_manual(values=c('white', 'black')) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position='none')
}
