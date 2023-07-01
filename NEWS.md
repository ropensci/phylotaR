# phylotaR 1.3.0

* Included code of `treeman` package in this package.
* Wrapped `outsider` package as it has been archived.
* Updated and fixed the dependency to `restez`.

There are a lot of CRAN submission problems pointed by Victoria Wimmer and I have no patience to handle them, I am not able to release it to CRAN.

```
Please add \value to .Rd files regarding exported methods and explain 
the functions results in the documentation. Please write about the 
structure of the output (class) and also what the output means. (If a 
function does not return a value, please document that too, e.g. 
\value{No return value, called for side effects} or similar)
Missing Rd-tags in up to 93 .Rd files, e.g.:
      addClade.Rd: \value
      addNdmtrx.Rd: \value
      addTip.Rd: \value
      blncdTree.Rd: \value
      calcDstBLD.Rd: \value
      calcDstMtrx.Rd: \value
      ...

You write information messages to the console that cannot be easily 
suppressed. It is more R like to generate objects that can be used to 
extract the information a user is interested in, and then print() that 
object.
Instead of print()/cat() rather use message()/warning()  or 
if(verbose)cat(..) (or maybe stop()) if you really have to write text to 
the console.
(except for print, summary, interactive functions)

Please do not modify the global environment (e.g. by using <<-) in your 
functions. This is not allowed by the CRAN policies.

Please fix and resubmit.
```

# phylotaR 1.2.0

## `outsider` integration

* Now work with `outsider`, reduing user need to install BLAST+
* New parameter: `db_only`, avoid downloads by fetching exclusively from
`restez` database
* Excluding Transcriptome Shotgun Assembly sequences, making downloading faster

# phylotaR 1.1.1

## Minor fixes

* Validity check for `SeqRec`

# phylotaR 1.1.0

## restez integration

* Now works with [`restez`](https://ropensci.github.io/restez/articles/4_phylotar.html)
* Multiple IDs in `txid` argument
* Fancy new colours
* Easier process killing

# phylotaR 1.0.0
