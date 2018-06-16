# Contributing

You are very welcome to help out in the development of phylotaR. If you have any ideas for future features than please add them to the [issues page](https://github.com/AntonelliLab/phylotaR/issues). If you have the impetus to add those features yourself, then please fork and send a pull request.

## Areas for possible contribution

### Alternatives to BLAST

Currently, phylotaR only makes use of BLAST. BLAST is good because it can work with any sequence type and is very sensitive. It is, however, slower than alternatives local alignment search tools. A set of generic input, run and output functions for working with any BLAST-alternative with an accompanying vignette would allow a user to any of the available alternatives. 

### BLAST API

A big issue with phylotaR is the need to install and run a local copy of BLAST. One possibility is add API functionality so that a user can send the jobs to the cloud instead of having to install their own version of BLAST.

### Inputting one's own sequences

A user may wish to make use of sequences they have generated themselves in conjunction with those available from GenBank. Currently there is no effort to allow a user to do this. It is a little tricky to do because of phylotaR's reliance on IDs.

### Multiple taxids

Many users wish to run a single phylotaR run for multiple, potentially paraphyletic, taxonomic IDs.

### A user determined taxonomy

Currrently, phylotaR depends on NCBI's taxonomy. In theory, a user could provide their own newick tree representing their preferred taxonomy instead.

### RefSeq

Identify sequences that are orthologous to a specified sequence.

## How to contribute

To contribute you will need a GitHub account and to have basic knowledge of the R language. You can then create a fork of the
repo in your own GitHub account and download the repository to your local machine. `devtools` is recommended.

```r
devtools::install_github('[your account]/phylotaR')
```

All new functions must be tested. For every new file in `R/`, a new test file must be created in `tests/testthat/`. To test the
package and make sure it meets CRAN guidelines use `devtools`. 

```r
devtools::test()
devtools::check_cran()
```

For help, refer to Hadley Wickham's book, [R packages](http://r-pkgs.had.co.nz/).

## Style guide

phylotaR is being developed for submission to ROpenSci. This means the package and its code should meet ROpenSci style and
standards. For example, function names should be all lowercase, separated by underscores and the last word should, ideally, be
a verb.

```
# e.g.
species_ids_retrieve()  # good
sppIDs()                # not great
sp.IDS_part2()          # really bad
sigNXTprt.p()           # awful
```

It is best to make functions small, with specific names. Feel free to break up code into multiple separate files (e.g. tools,
helper functions, stages ...). For more details and better explanations refer to the ROpenSci [guide](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md).
