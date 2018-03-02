#' @name calcClstrs
#' @title Calculate clusters for all sequences in WD
#' @description TODO
calcClstrs <- function(txdct, ps) {
  # load sequences
  sq_fls <- list.files(file.path(ps[['wd']], 'cache', 'sqs'))
  #foreach(i=seq_along(sq_fls)) %dopar% {
  for(i in seq_along(sq_fls)) {
    sq_fl <- sq_fls[i]
    # TODO: use the cache tool
    sqs <- readRDS(file=file.path(file.path(ps[['wd']], 'cache',
                                            'sqs', sq_fl)))
    txid <- as.character(sub('\\.RData', '', sq_fl))
    info(lvl=1, ps=ps, "Working on [id ", txid, "]")
    clstrs <- clstrAll(txid=txid, sqs=sqs, txdct=txdct,
                       ps=ps)
    svClstrs(wd=ps[['wd']], txid=txid, clstrs=clstrs)
    info(lvl=1, ps=ps, "Finished [id ", txid, "] # [",
        i, "/", length(sq_fls), "]")
  }
}
