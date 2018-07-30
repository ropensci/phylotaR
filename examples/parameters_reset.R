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
  # run
  # run(wd = wd) # not running in test
  # use ctrl+c or Esc to kill
  # change parameters, e.g. min and max sequence lengths
  parameters_reset(wd = 'aotus', parameters = c('mnsql', 'mxsql'),
                   values = c(300, 1500))
  # see ?parameters
  # restart
  restart(wd = wd)
}
