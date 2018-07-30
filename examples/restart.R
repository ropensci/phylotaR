\dontrun{
  
  # Note: this example requires BLAST and internet to run.
  
  # example with temp folder
  wd <- file.path(tempdir(), 'aotus')
  # setup for aotus, make sure aotus/ folder already exists
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  ncbi_dr <- '[SET BLAST+ BIN PATH HERE]'
  setup(wd = wd, txid = 9504, ncbi_dr = ncbi_dr)  # txid for Aotus primate genus
  # run and stop after 10 seconds
  R.utils::withTimeout(expr = {
    run(wd = wd)
  }, timeout = 10)
  # use ctrl+c or Esc to kill without a timelimit
  # and restart with ....
  restart(wd = wd)
}
