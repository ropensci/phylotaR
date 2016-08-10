## script to populate the 'accession2taxid' table, needed
## to create the other phylota tables

accession2taxid.create <- function(dbloc, taxdir) {
    subdir <- 'accession2taxid'
    file.loc <- paste0(taxdir, '/', subdir, '/nucl_gb.accession2taxid')    
    cmd <- paste('sqlite3', dbloc, '\'.import', file.loc, 'accession2taxid\'')
    system(cmd)
    ## create index    
    icmd <- '\'create index taxid on accession2taxid(taxid)\'';
}
