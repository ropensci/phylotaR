\dontrun{
  
  # reset and restart aotus run
  setup(wd = 'aotus', txid = 9504)
  run(wd = 'aotus')
  # use ctrl+c or Esc to kill
  # reset to an earlier stage
  reset(wd = 'aotus', stage = 'download')
  # restart with from download
  restart(wd = 'aotus')
}
