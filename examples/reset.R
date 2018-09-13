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
  # individually run taxise
  taxise_run(wd = wd)
  # reset back to taxise as if it has not been run
  reset(wd = 'aotus', stage = 'taxise')
  # run taxise again ....
  taxise_run(wd = wd)
}
