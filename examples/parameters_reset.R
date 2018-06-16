\dontrun{
  
  # reset and restart aotus run
  setup(wd = 'aotus', txid = 9504)
  run(wd = 'aotus')
  # use ctrl+c or Esc to kill
  # change parameters, e.g. reduce btchsz to 100 if NCBI is blocking
  parameters_reset(wd = 'aotus', parameters = 'btchsz', values = 100)
  # restart
  restart(wd = 'aotus')
}
