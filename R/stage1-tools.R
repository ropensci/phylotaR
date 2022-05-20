#' @name tax_download
#' @title Download taxonomic records
#' @description Downloads one batch of taxonomic
#' records.
#' @return list of list
#' @param ids Vector of taxonomic IDs
#' @template ps
#' @family run-private
tax_download <- function(ids, ps) {
  recs <- vector("list", length = length(ids))
  names(recs) <- ids
  args <- list(db = "taxonomy", id = ids, rettype = "xml")
  raw_recs <- search_and_cache(
    func = rentrez::entrez_fetch, fnm = "fetch",
    args = args, ps = ps
  )
  raw_recs <- XML::xmlToList(raw_recs)
  if (length(raw_recs) != length(ids)) {
    # remove cached object
    ncbicache_save(fnm = "fetch", args = args, wd = ps[["wd"]], obj = NULL)
    warn(
      ps, "Not all taxonomic records could be downloaded.",
      "NCBI may be limiting usage, run taxise step again."
    )
    return(NULL)
  }
  for (i in seq_along(raw_recs)) {
    raw_rec <- raw_recs[[i]]
    rw_lng <- raw_rec[["LineageEx"]]
    lng_ids <- vapply(rw_lng, function(x) x[["TaxId"]], "")
    names(lng_ids) <- NULL
    lng_ids <- c(lng_ids, raw_rec[["TaxId"]])
    lng_rnks <- vapply(rw_lng, function(x) x[["Rank"]], "")
    names(lng_rnks) <- NULL
    lng_rnks <- c(lng_rnks, raw_rec[["Rank"]])
    lng <- list("ids" = lng_ids, "rnks" = lng_rnks)
    cmnms <- raw_rec[["OtherNames"]]
    if ("GenbankCommonName" %in% names(cmnms)) {
      cmnm <- cmnms[["GenbankCommonName"]]
    } else if ("CommonName" %in% names(cmnms)) {
      cmnm <- cmnms[["CommonName"]]
    } else {
      cmnm <- ""
    }
    rec <- new("TaxRec",
      id = raw_rec[["TaxId"]],
      scnm = raw_rec[["ScientificName"]],
      cmnm = cmnm, rnk = raw_rec[["Rank"]],
      lng = lng, prnt = raw_rec[["ParentTaxId"]]
    )
    recs[[rec@id]] <- rec
  }
  recs
}

parent_recs_get <- function(recs, txids, ps) {
  # look up the root taxon for all txids provided
  lngs <- lapply(X = recs[ps[["txid"]]], FUN = function(x) {
    slot(object = x, name = "lng")[["ids"]]
  })
  candidates <- lngs[[1]]
  for (lng in lngs) {
    candidates <- candidates[candidates %in% lng]
  }
  root_txid <- candidates[length(candidates)]
  node_txids <- unique(unname(unlist(lapply(X = lngs, function(x) {
    x[!x %in% candidates]
  }))))
  node_txids <- node_txids[!node_txids %in% txids]
  node_txids <- c(node_txids, root_txid)
  obj_save(wd = ps[["wd"]], obj = node_txids, nm = "ignore_ids")
  node_recs <- batcher(ids = node_txids, func = tax_download, ps = ps, lvl = 2)
  recs <- c(recs, node_recs)
  txids <- c(txids, node_txids)
  list("recs" = recs, "txids" = txids)
}
