# treeman package -------------------------------
#' @name calcNdBlnc
#' @title Calculate the balance of a node
#' @description Returns the balance of a node.
#' @details Balance is calculated as the absolute difference between the number of descendents
#' of the two bifurcating edges of a node and the expected value for a balanced tree.
#' \code{NA} is returned if the node is polytomous or a tip.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{calcNdsBlnc}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' calcNdBlnc(tree, id = tree["root"]) # root balance
calcNdBlnc <- function(tree, id) {
  ntot <- length(getNdKids(tree, id))
  ptids <- tree@ndlst[[id]][["ptid"]]
  if (length(ptids) > 2) {
    return(NA)
  }
  ptid <- ptids[1]
  nprt <- length(getNdKids(tree, ptid))
  if (nprt == 0) {
    nprt <- 1
  }
  abs((ntot / 2) - nprt)
}

#' @name calcNdsBlnc
#' @title Calculate the balances of all nodes
#' @description Returns the absolute differences in number of descendants for bifurcating
#' branches of every node
#' @details Runs \code{calcNdBlnc()} across all node IDs. \code{NA} is returned if the
#' node is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{calcNdBlnc}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' calcNdsBlnc(tree, ids = tree["nds"])
calcNdsBlnc <- function(tree, ids, parallel = FALSE, progress = "none") {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  plyr::mdply(
    .data = l_data, .fun = calcNdBlnc, tree = tree,
    .parallel = parallel, .progress = progress
  )[, 2]
}

#' @name calcDstTrp
#' @title Calculate the triplet distance between two trees
#' @description Returns the triplet distance between two trees.
#' @details The triplet distance is calculated as the sum of different outgroups among
#' every triplet of tips between the two trees. Normalisation is performed by dividing
#' the resulting number by the total number of triplets shared between the two trees.
#' The triplet distance is calculated only for shared tips between the two trees. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Critchlow DE, Pearl DK, Qian C. (1996) The Triples Distance for rooted bifurcating phylogenetic trees.
#' Systematic Biology, 45, 323-34.
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstRF}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstTrp(tree_1, tree_2)
calcDstTrp <- function(tree_1, tree_2, nrmlsd = FALSE,
                       parallel = FALSE, progress = "none") {
  .count <- function(i) {
    o1 <- getOtgrp(tree_1, cmbs[, i])
    o2 <- getOtgrp(tree_2, cmbs[, i])
    if (length(o1) != length(o2) || o1 != o2) {
      cntr <- 1
    } else {
      cntr <- 0
    }
    cntr
  }
  shrd <- tree_1@tips[tree_1@tips %in% tree_2@tips]
  cmbs <- combn(shrd, 3)
  l_data <- data.frame(i = 1:ncol(cmbs), stringsAsFactors = FALSE)
  res <- plyr::mdply(.data = l_data, .count, .parallel = parallel, .progress = progress)
  cntr <- sum(res[, 2])
  if (nrmlsd) {
    cntr <- cntr / ncol(cmbs)
  }
  cntr
}

#' @name calcOvrlp
#' @title Calculate phylogenetic overlap
#' @description Returns the sum of branch lengths represented by ids_1 and ids_2 for a tree.
#' @details Use this to calculate the sum of branch lengths that are represented between two
#' communities. This measure is also known as the unique fraction. It can be used to measure
#' concepts of phylogenetic turnover. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids_1 tip ids of community 1
#' @param ids_2 tip ids of community 2
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing
#' microbial communities. Applied and Environmental Microbiology, 71(12), 8228-35.
#' @seealso
#' \code{\link{calcPhyDv}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' ids_1 <- sample(tree["tips"], 5)
#' ids_2 <- sample(tree["tips"], 5)
#' calcOvrlp(tree, ids_1, ids_2)
calcOvrlp <- function(tree, ids_1, ids_2, nrmlsd = FALSE,
                      parallel = FALSE, progress = "none") {
  if (progress != "none") {
    cat("Part 1/3 ....\n")
  }
  spns <- getNdsSlt(tree,
    slt_nm = "spn", ids = tree@all,
    parallel = parallel, progress = progress
  )
  names(spns) <- tree@all
  if (progress != "none") {
    cat("Part 2/3 ....\n")
  }
  ids_1 <- c(unique(unlist(getNdsPrids(tree,
    ids = ids_1,
    parallel = parallel,
    progress = progress
  ))), ids_1)
  if (progress != "none") {
    cat("Part 3/3 ....\n")
  }
  ids_2 <- c(unique(unlist(getNdsPrids(tree,
    ids = ids_2,
    parallel = parallel,
    progress = progress
  ))), ids_2)
  ovrlp <- sum(spns[ids_2[ids_2 %in% ids_1]])
  if (nrmlsd) {
    ovrlp <- ovrlp / tree@pd
  }
  ovrlp
}

#' @name calcDstBLD
#' @title Calculate the BLD between two trees
#' @description Returns the branch length distance between two trees.
#' @details BLD is the Robinson-Foulds distance weighted by branch length. Instead of summing
#' the differences in partitions between the two trees, the metric takes the square root
#' of the squared difference in branch lengths. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Kuhner, M. K. and Felsenstein, J. (1994) Simulation comparison of phylogeny
#' algorithms under equal and unequal evolutionary rates. Molecular Biology and
#' Evolution, 11, 459-468.
#' @seealso
#' \code{\link{calcDstTrp}}, \code{\link{calcDstRF}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstBLD(tree_1, tree_2)
calcDstBLD <- function(tree_1, tree_2, nrmlsd = FALSE,
                       parallel = FALSE, progress = "none") {
  n1 <- tree_1@nds[!tree_1@nds == tree_1@root]
  n2 <- tree_2@nds[!tree_2@nds == tree_2@root]
  if (progress != "none") {
    cat("Part 1/2 ....\n")
  }
  c1 <- getNdsKids(tree_1, n1, parallel = parallel, progress = progress)
  if (progress != "none") {
    cat("Part 2/2 ....\n")
  }
  c2 <- getNdsKids(tree_2, n2, parallel = parallel, progress = progress)
  s1 <- getNdsSlt(tree_1, slt_nm = "spn", ids = n1)
  s2 <- getNdsSlt(tree_2, slt_nm = "spn", ids = n2)
  d1 <- s2[match(c1, c2)]
  d1[which(is.na(d1))] <- 0
  d1 <- s1 - d1
  d2 <- s1[match(c2, c1)]
  d2[which(is.na(d2))] <- 0
  d2 <- s2 - d2
  d <- sqrt(sum(c(d1^2, d2^2)))
  if (nrmlsd) {
    max_d <- sqrt(sum(c(s1^2, s2^2)))
    d <- d / max_d
  }
  d
}

#' @name calcDstRF
#' @title Calculate the Robinson-Foulds distance between two trees
#' @description Returns the Robinson-Foulds distance between two trees.
#' @details RF distance is calculated as the sum of partitions in one tree that are
#' not shared by the other. The maximum number of split differences is the total number
#' of nodes in both trees (excluding the roots). Trees are assumed to be bifurcating,
#' this is not tested. The metric is calculated as if trees are unrooted. Parallelizable.
#' @param tree_1 \code{TreeMan} object
#' @param tree_2 \code{TreeMan} object
#' @param nrmlsd Boolean, should returned value be between 0 and 1? Default, FALSE.
#' @references
#' Robinson, D. R.; Foulds, L. R. (1981). "Comparison of phylogenetic trees".
#' Mathematical Biosciences 53: 131-147.
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstTrp}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' calcDstRF(tree_1, tree_2)
calcDstRF <- function(tree_1, tree_2, nrmlsd = FALSE) {
  shared_tips <- c(tree_1@tips, tree_2@tips)
  shared_tips <- shared_tips[duplicated(shared_tips)]
  b1 <- getBiprts(
    tree = tree_1, tips = shared_tips, universal = TRUE,
    root = FALSE
  )
  b2 <- getBiprts(
    tree = tree_2, tips = shared_tips, universal = TRUE,
    root = FALSE
  )
  # count unique paritions and sum to calc RFD
  d <- sum(!b1 %in% b2) + sum(!b2 %in% b1)
  if (nrmlsd) {
    max_d <- (length(b1) + length(b2))
    d <- d / max_d
  }
  d
}

#' @name calcPhyDv
#' @title Calculate phylogenetic diversity
#' @description Returns the phylogenetic diversity of a tree for the tips specified.
#' @details Faith's phylogenetic diversity is calculated as the sum of all connected
#' branches for specified tips in a tree. It can be used to investigate how biodviersity
#' as measured by the phylogeny changes. Parallelizable.
#' The function uses \code{getCnntdNds()}.
#' @param tree \code{TreeMan} object
#' @param tids tip ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Faith, D. (1992). Conservation evaluation and phylogenetic diversity.
#'  Biological Conservation, 61, 1-10.
#' @seealso
#' \code{\link{calcFrPrp}}, \code{\link{calcOvrlp}}, \code{\link{getCnnctdNds}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' calcPhyDv(tree, tree["tips"])
calcPhyDv <- function(tree, tids,
                      parallel = FALSE, progress = "none") {
  prids <- getCnnctdNds(tree, tids)
  spns <- getNdsSlt(tree,
    slt_nm = "spn", ids = prids,
    parallel = parallel, progress = progress
  )
  sum(spns)
}

#' @name calcFrPrp
#' @title Calculate evolutionary distinctness
#' @description Returns the evolutationary distinctness of ids using the fair proportion metric.
#' @details The fair proportion metric calculates the evolutionary distinctness of tips
#' in a tree through summing the total amount of branch length each tip represents, where
#' each branch in the tree is evenly divided between all descendants. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007).
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#' @seealso
#' \code{\link{calcPhyDv}}, \code{\link{calcPrtFrPrp}},
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' calcFrPrp(tree, tree["tips"])
calcFrPrp <- function(tree, tids, progress = "none") {
  .calc <- function(i) {
    id <- tree@all[i]
    spn <- getNdSlt(tree, "spn", id)
    if (id %in% tree@tips) {
      spn_shres[i, id] <<- spn
    } else {
      kids <- getNdKids(tree, id)
      spn_shre <- spn / length(kids)
      spn_shres[i, kids] <<- spn_shre
    }
  }
  spn_shres <- matrix(0, ncol = tree@ntips, nrow = tree@nall)
  colnames(spn_shres) <- tree@tips
  plyr::m_ply(
    .data = data.frame(i = 1:tree@nall), .fun = .calc,
    .progress = progress
  )
  colSums(spn_shres[, tids])
}

#' @name calcPrtFrPrp
#' @title Calculate evolutionary distinctness for part of tree
#' @description Returns the evolutationary distinctness of ids using the fair proportion metric.
#' @details Extension of \code{calcFrPrp()} but with ignore argument.
#' Use \code{ignr} to ignore certain tips from calculation. For example, if any of tips
#' are extinct you may wish to ignore these.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param ignr tips to ignore in calculation
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @references
#' Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007).
#'  Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#' @seealso
#' \code{\link{calcFrPrp}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' calcPrtFrPrp(tree, c("t1", "t3"), ignr = "t2")
calcPrtFrPrp <- function(tree, tids, ignr = NULL, progress = "none") {
  .calc <- function(i) {
    id <- allnds[i]
    spn <- getNdSlt(tree, "spn", id)
    if (id %in% tips) {
      spn_shres[i, id] <<- spn
    } else {
      kids <- getNdKids(tree, id)
      kids <- kids[!kids %in% ignr]
      if (length(kids) > 0) {
        spn_shre <- spn / length(kids)
        spn_shres[i, kids] <<- spn_shre
      }
    }
  }
  tips <- tree@tips[!tree@tips %in% ignr]
  allnds <- tree@all[!tree@all %in% ignr]
  spn_shres <- matrix(0, ncol = length(tips), nrow = length(allnds))
  colnames(spn_shres) <- tips
  plyr::m_ply(
    .data = data.frame(i = 1:length(allnds)), .fun = .calc,
    .progress = progress
  )
  colSums(spn_shres[, tids])
}

#' @name calcDstMtrx
#' @title Calculate the distance matrix
#' @description Returns a distance matrix for specified ids of a tree.
#' @details The distance between every id in the tree is calculated by summing the
#' lengths of the branches that connect them. This can be useful for testing the distances
#' between trees, checking for evoltuionary isolated tips etc. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids IDs of nodes/tips
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{calcDstBLD}}, \code{\link{calcDstRF}}, \code{\link{calcDstTrp}}
#' \url{https://github.com/DomBennett/treeman/wiki/calc-methods}
#' @export
#' @examples
#' # checking the distance between two trees
#'
#' tree_1 <- randTree(10)
#' tree_2 <- randTree(10)
#' dmat1 <- calcDstMtrx(tree_1, tree_1["tips"])
#' dmat2 <- calcDstMtrx(tree_2, tree_2["tips"])
#' mdl <- cor.test(x = dmat1, y = dmat2)
#' as.numeric(1 - mdl$estimate) # 1 - Pearson's r
calcDstMtrx <- function(tree, ids, parallel = FALSE,
                        progress = "none") {
  .getDist <- function(id_1, id_2) {
    if (id_1 == id_2) {
      return(0)
    }
    path <- getPath(tree, from = id_1, to = id_2)
    path_spns <- unlist(lapply(tree@ndlst[path], function(n) n[["spn"]]))
    sum(path_spns)
  }
  cmbs <- expand.grid(ids, ids, stringsAsFactors = FALSE)
  colnames(cmbs) <- c("id_1", "id_2")
  res <- plyr::mdply(
    .data = cmbs, .fun = .getDist,
    .parallel = parallel, .progress = progress
  )
  res <- matrix(res[, 3], ncol = length(ids))
  colnames(res) <- rownames(res) <- ids
  res
}

#' @name fastCheckTreeMan
#' @title Check if tree is correct, fast!
#' @description Return T/F if tree is a true \code{TreeMan} object
#' @details Whenever a tree is first initiated this check is used.
#' For more detailed checking use \code{checkNdlst}.
#' @param object \code{TreeMan} object
#' @seealso
#' \code{\link{checkNdlst}}, \code{\link{checkTreeMen}}
#' @export
fastCheckTreeMan <- function(object) {
  kwn_ids <- names(object@ndlst)
  ids <- unlist(sapply(object@ndlst, function(x) x[c("id", "ptid", "prid")]))
  test <- all(ids %in% kwn_ids) & object@root %in% kwn_ids
  # check hierarchy through prinds
  prinds <- object@prinds
  if (length(prinds) > 0) {
    # only root node should equal its index
    prind_test <- sum(prinds == 1:length(prinds)) == 1
    # all internal nodes should occur more than once (twice for bifurcating trees)
    prind_test <- all(table(prinds) > 1) & prind_test
    test <- test & prind_test
  }
  test
}

#' @name checkNdlst
#' @title Check if ndlst is correct
#' @description Return T/F fpr \code{ndlst} consistency
#' @details Tests whether each node in tree points to valid other node IDs. Also
#' ensures `spn` and `root` are correct. Reports nodes that have errors.
#' @param ndlst \code{ndlst}
#' @param root root ID
#' @seealso
#' \code{\link{fastCheckTreeMan}}, \code{\link{checkTreeMen}}
#' @export
#' @examples
#'
#' tree <- randTree(100)
#' (checkNdlst(tree@ndlst, tree@root))
checkNdlst <- function(ndlst, root) {
  .check <- function(nd) {
    # must have id
    test_id <- is.character(nd[["id"]]) & "id" %in% names(nd)
    # id must contain no special characters
    test_spcl_chrs <- test_id && !grepl("[^0-9a-zA-Z_]", nd[["id"]])
    # txnyms
    test_txnym <- TRUE
    if ("txnym" %in% names(nd)) {
      test_txnym <- is.character(nd[["txnym"]])
      for (txnym in nd[["txnym"]]) {
        test_txnym <- test_txnym && !grepl("[^0-9a-zA-Z_]", txnym)
      }
    }
    # must have either prid/ptid or both
    test_slts <- ("ptid" %in% names(nd) | "prid" %in% names(nd))
    test_valid_nd <- nd[["id"]] %in% nds # nd id must be known
    # prid and ptids must be known
    test_prid <- is.character(nd[["prid"]]) & nd[["prid"]] %in% nds
    test_ptid <- is.character(nd[["ptid"]]) > 0 & all(nd[["ptid"]] %in% nds)
    # spns must be 0 or more
    test_spn <- is.numeric(nd[["spn"]]) && nd[["spn"]] >= 0
    # test self-reference
    test_sr <- !nd[["prid"]] %in% nd[["ptid"]]
    # test root is never a ptid, proxy for circularity
    test_circ <- !rid %in% nd[["ptid"]]
    # only root is self-referential
    test_root <- rid != nd[["id"]] |
      (rid == nd[["id"]] & rid == nd[["prid"]])
    bool <- test_id & test_valid_nd & test_prid &
      test_ptid & test_sr & test_circ & test_slts &
      test_spcl_chrs & test_txnym
    if (length(bool) > 0 && bool) {
      return(TRUE)
    }
    FALSE
  }
  nds <- names(ndlst)
  if (any(duplicated(nds))) {
    dups <- nds[duplicated(nds)]
    dups <- unique(dups)
    msg <- "These node IDs are duplicated:\n"
    if (length(dups) > 1) {
      for (i in 1:length(dups) - 1) {
        msg <- paste0(msg, dups[i], ", ")
      }
    }
    msg <- paste0(msg, dups[length(dups)], "\n")
    cat(msg)
    return(FALSE)
  }
  rid <- root
  nd_checks <- vapply(ndlst, .check, logical(1))
  if (!all(nd_checks)) {
    msg <- "These nodes are invalid:\n"
    bad <- which(!nd_checks)
    if (length(bad) > 1) {
      for (i in bad[-length(bad)]) {
        msg <- paste0(msg, nds[i], ", ")
      }
    }
    msg <- paste0(msg, nds[bad[length(bad)]], "\n")
    cat(msg, "\n")
    return(FALSE)
  }
  if (!rid %in% nds) {
    msg <- paste0("Root node `", rid, "` not in ndlst\n")
    cat(msg, "\n")
    return(FALSE)
  }
  TRUE
}

#' @name checkTreeMen
#' @title Check if trees are correct
#' @description Return T/F if trees is a true \code{TreeMen} object
#' @details Tests whether all trees in object are \code{TreeMan} objects
#' @param object \code{TreeMen} object
#' @seealso
#' \code{\link{checkNdlst}}
#' @export
checkTreeMen <- function(object) {
  .check <- function(i) {
    if (class(object@treelst[[i]])[1] != "TreeMan") {
      invlds <<- c(i, invlds)
    }
    NULL
  }
  invlds <- NULL
  mapply(.check, i = 1:length(object@treelst))
  if (length(invlds) > 0) {
    for (i in invlds) {
      cat("[", i, "] in treelst is not a TreeMan object\n", sep = "")
    }
    return(FALSE)
  }
  TRUE
}

#' phylo class
#'
#' @name phylo-class
#' @aliases phylo
#'
#' @exportClass phylo
setOldClass("phylo")

#' multiPhylo class
#'
#' @name multiPhylo-class
#' @aliases multiPhylo
#'
#' @exportClass multiPhylo
setOldClass("multiPhylo")



#' @name randTree
#' @title Generate a random tree
#' @description Returns a random \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{rtree()} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{blncdTree}},
#' \code{\link{unblncdTree}}
#' @export
#' @examples
#'
#' tree <- randTree(5)
randTree <- function(n, wndmtrx = FALSE, parallel = FALSE) {
  # Return a random tree based on a broken-stick model
  .randomPrinds <- function(n) {
    pool <- rep((1:(n - 1)), each = 2)
    res <- rep(NA, length(pool) + 1)
    res[1] <- 1
    for (i in 2:length(res)) {
      pssbls <- which(i > pool)
      if (length(pssbls) == 1) {
        i_pool <- pssbls
      } else {
        i_pool <- sample(pssbls, 1)
      }
      res[i] <- pool[i_pool]
      pool[i_pool] <- NA
    }
    res
  }
  if (n < 3) {
    stop("`n` is too small")
  }
  prinds <- .randomPrinds(n)
  .cnstrctTree(n, prinds,
    wndmtrx = wndmtrx,
    parallel = parallel
  )
}

#' @name blncdTree
#' @title Generate a balanced tree
#' @description Returns a balanced \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{stree(type='balanced')} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}},
#' \code{\link{unblncdTree}}
#' @export
#' @examples
#'
#' tree <- blncdTree(5)
blncdTree <- function(n, wndmtrx = FALSE, parallel = FALSE) {
  if (n < 3) {
    stop("`n` is too small")
  }
  prinds <- c(1, rep((1:(n - 1)), each = 2))
  .cnstrctTree(n, prinds,
    wndmtrx = wndmtrx,
    parallel = parallel
  )
}

#' @name unblncdTree
#' @title Generate an unbalanced tree
#' @description Returns an unbalanced \code{TreeMan} tree with \code{n}
#' tips.
#' @details Equivalent to \code{ape}'s \code{stree(type='left')} but returns a
#' \code{TreeMan} tree. Tree is always rooted and bifurcating.
#' @param n number of tips, integer, must be 3 or greater
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel T/F run in parallel? Default FALSE.
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}},
#' \code{\link{blncdTree}}
#' @export
#' @examples
#'
#' tree <- unblncdTree(5)
unblncdTree <- function(n, wndmtrx = FALSE, parallel = FALSE) {
  if (n < 3) {
    stop("`n` is too small")
  }
  prinds <- c(1, 1:(n - 1), 1:(n - 1))
  .cnstrctTree(n, prinds,
    wndmtrx = wndmtrx,
    parallel = parallel
  )
}

.cnstrctTree <- function(n, prinds, wndmtrx, parallel) {
  .add <- function(i) {
    nd <- vector("list", length = 4)
    names(nd) <- c("id", "ptid", "prid", "spn")
    nd[["id"]] <- ids[i]
    nd[["spn"]] <- spns[i]
    nd[["prid"]] <- ids[prinds[i]]
    nd[["ptid"]] <- ptids[ptnds_pool == i]
    nd
  }
  nnds <- length(prinds)
  # random numbers for spans
  spns <- c(0, runif(nnds - 1, 0, 1))
  ids <- rep(NA, nnds)
  tinds <- which(!1:nnds %in% prinds)
  ids[tinds] <- paste0("t", 1:n)
  ids[1:nnds %in% prinds] <- paste0("n", 1:(n - 1))
  ptnds_pool <- prinds[-1]
  ptids <- ids[-1]
  ndlst <- plyr::mlply(.data = 1:nnds, .fun = .add, .parallel = parallel)
  attr(ndlst, "split_labels") <-
    attr(ndlst, "split_type") <- NULL
  names(ndlst) <- ids
  tree <- new("TreeMan",
    ndlst = ndlst, root = "n1", wtxnyms = FALSE,
    ndmtrx = NULL, prinds = prinds, tinds = tinds
  )
  tree <- updateSlts(tree)
  if (wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name twoer
#' @title Generate a tree of two tips
#' @description Returns a \code{TreeMan} tree with two tips and a root.
#' @details Useful for building larger trees with \code{addClade()}.
#' Note, a node matrix cannot be added to a tree of two tips.
#' @param tids tip IDs
#' @param spns tip spans
#' @param rid root ID
#' @param root_spn root span
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{randTree}}
#' @export
#' @examples
#'
#' tree <- twoer()
twoer <- function(tids = c("t1", "t2"), spns = c(1, 1),
                  rid = "root", root_spn = 0) {
  ndlst <- list()
  ndlst[[rid]] <- list(
    "id" = rid, "prid" = rid,
    "ptid" = tids[1:2], "spn" = root_spn
  )
  ndlst[[tids[1]]] <- list(
    "id" = tids[[1]], "prid" = rid,
    "ptid" = NULL, "spn" = spns[1]
  )
  ndlst[[tids[2]]] <- list(
    "id" = tids[[2]], "prid" = rid,
    "ptid" = NULL, "spn" = spns[2]
  )
  prinds <- c(1, 1, 1)
  tinds <- c(2, 3)
  tree <- new("TreeMan",
    ndlst = ndlst, root = "root", wtxnyms = FALSE,
    ndmtrx = NULL, prinds = prinds, tinds = tinds
  )
  updateSlts(tree)
}

#' @name getNdLng
#' @title Get lineage
#' @description Return unique taxonomic names for connecting \code{id} to root.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsLng}}, \code{\link{getNdsFrmTxnyms}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # return human lineage
#' getNdLng(mammals, id = "Homo_sapiens")
getNdLng <- function(tree, id) {
  .get <- function(txnym, ...) {
    lng <<- c(txnym, lng)
  }
  prids <- getNdPrids(tree, id)
  lng <- NULL
  plyr::m_ply(tree@ndlst[prids], .fun = .get)
  unique(lng)
}

#' @name getNdAge
#' @title Get age
#' @description Return the age for \code{id}. Requires the known age of the tree to be provided.
#' @details Returns a numeric.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree
#' @seealso
#' \code{\link{getNdsAge}},
#' \code{\link{getSpnAge}},
#' \code{\link{getSpnsAge}},
#' \code{\link{getPrnt}}, \code{\link{getAge}}
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # when did apes emerge?
#' # get parent id for all apes
#' prnt_id <- getPrnt(mammals, ids = c("Homo_sapiens", "Hylobates_concolor"))
#' # mammal_age <- getAge(mammals)  # ~166.2, needs to be performed when tree is not up-to-date
#' getNdAge(mammals, id = prnt_id, tree_age = 166.2)
getNdAge <- function(tree, id, tree_age) {
  tree_age - .getNdPrdstsFrmLst(tree@ndlst, tree@prinds, id = id)
}

#' @name getSpnAge
#' @title Get age range
#' @description Return start and end ages for \code{id} from when it first appears to when it splits
#' @details Returns a dataframe.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @param tree_age numeric value of known age of tree
#' @seealso
#' \code{\link{getNdAge}},
#' \code{\link{getNdsAge}},
#' \code{\link{getSpnsAge}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # mammal_age <- getAge(mammals)  # ~166.2, needs to be performed when tree is not up-to-date
#' getSpnAge(mammals, id = "Homo_sapiens", tree_age = 166.2)
getSpnAge <- function(tree, id, tree_age) {
  end <- .getNdPrdstsFrmLst(tree@ndlst, tree@prinds, id = id)
  start <- end - tree@ndlst[[id]][["spn"]]
  end <- tree_age - end
  start <- tree_age - start
  data.frame(spn = id, start, end)
}

#' @name getNdPrids
#' @title Get pre-nodes to root
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a vector. IDs are returned order from node ID to root.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrids}},
#' \code{\link{getNdPtids}},
#' \code{\link{getNdsPtids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # get all nodes to root
#' getNdPrids(tree, id = "t1")
getNdPrids <- function(tree, id) {
  .getNdPridsFrmLst(tree@ndlst, prinds = tree@prinds, id = id)
}

#' @name getNdPtids
#' @title Get post-nodes to tips
#' @description Return node ids for connecting \code{id} to kids.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPtids}},
#' \code{\link{getNdPrids}},
#' \code{\link{getNdsPrids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # get all nodes from root to tip
#' getNdPtids(tree, id = "n1")
# reduce dependence on the recursive, by getting prenodes
# tip ids to id
getNdPtids <- function(tree, id) {
  .getNdPtidsFrmLst(tree@ndlst, tree@prinds, id = id)
}

#' @name getNdKids
#' @title Get children IDs
#' @description Return the node ids of all tips that descend from node.
#' @details Returns a vector
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsKids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # everyone descends from root
#' getNdKids(tree, id = tree["root"])
getNdKids <- function(tree, id) {
  .getNdKidsFrmLst(tree@ndlst, prinds = tree@prinds, tinds = tree@tinds, id = id)
}

#' @name getNdPrdst
#' @title Get pre-distance
#' @description Return root to tip distance (prdst) for \code{id}
#' @details Sums the lengths of all branches from \code{id} to root.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPrdst}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdPrdst(tree, id = "t1") # return the distance to root from t1
getNdPrdst <- function(tree, id) {
  .getNdPrdstsFrmLst(tree@ndlst, tree@prinds, id = id)
}

#' @name getNdSlt
#' @title Get a node slot
#' @description Returns the value of named slot.
#' @details Returned object depends on name, either character, vector or numeric.
#' Default node slots are: id, spn, prid, ptid and txnym. If slot is empty, returns NA.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param id node id
#' @seealso
#' \code{\link{getNdsSlt}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdSlt(tree, slt_nm = "spn", id = "t1") # return span of t1
getNdSlt <- function(tree, slt_nm, id) {
  res <- tree@ndlst[[id]][[slt_nm]]
  if (is.null(res)) {
    res <- NA
  }
  res
}

#' @name getNdPD
#' @title Get phylogenetic diversity of node
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsPD}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdPD(tree, id = "n1") # return PD of n1 which in this case is for the whole tree
getNdPD <- function(tree, id) {
  .getNdPDFrmLst(tree@ndlst, prinds = tree@prinds, id = id)
}

#' @name getNdSstr
#' @title Get sister id
#' @description Returns the id of the sister(s) of node id given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getNdsSstr}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdSstr(tree, id = "t1")
getNdSstr <- function(tree, id) {
  if (id == tree@root) {
    return(NULL)
  }
  .getNdSstrFrmLst(tree@ndlst, id)
}

#' @name getNdsLng
#' @title Get lineage for multiple nodes
#' @description Return unique taxonyms for connecting \code{ids} to root.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdLng}}, \code{\link{getNdsFrmTxnyms}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # return human and gorilla lineages
#' getNdsLng(mammals, id = c("Homo_sapiens", "Gorilla_gorilla"))
getNdsLng <- function(tree, ids, parallel = FALSE, progress = "none") {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  out <- plyr::mlply(
    .data = l_data, .fun = getNdLng, tree = tree,
    .parallel = parallel, .progress = progress
  )
  names(out) <- attr(out, "split_labels")[, 1]
  res <- out[1:length(out)]
  res
}

#' @name getNdsSstr
#' @title Get sister id
#' @description Returns the ids of the sister(s) of nd ids given.
#' @details An error is raised if there is no sister (e.g. for the root).
#'  There can be more than one sister if tree is polytomous. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids nd ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdSstr}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsSstr(tree, ids = tree["tips"])
getNdsSstr <- function(tree, ids, parallel = FALSE, progress = "none") {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  res <- plyr::mdply(
    .data = l_data, .fun = .getNdSstrFrmLst, ndlst = tree@ndlst,
    .parallel = parallel, .progress = progress
  )
  res[, 2]
}

#' @name getNdsPD
#' @title Get phylogenetic diversities of nodes
#' @description Return summed value of all descending spns
#' @details Sums the lengths of all descending branches from a node.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPD}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsPD(tree, ids = tree["all"]) # return PD of all ids
getNdsPD <- function(tree, ids, parallel = FALSE, progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    all_ids <- tree@all
    spns <- .getSltSpns(tree@ndlst)
    res <- .getNdsPDFrmMtrx(
      tree@ndmtrx, all_ids, ids, spns,
      parallel, progress
    )
  } else {
    res <- .getNdsPDFrmLst(tree@ndlst,
      prinds = tree@prinds,
      ids = ids, parallel = parallel, progress = progress
    )
  }
  res
}

#' @name getNdsPrdst
#' @title Get pre-distances
#' @description Return root to tip distances (prdst) for \code{ids}
#' @details Sums the lengths of all branches from \code{ids} to root.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPrdst}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsPrdst(tree, ids = tree["tips"]) # return prdsts for all tips
getNdsPrdst <- function(tree, ids, parallel = FALSE, progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    all_ids <- tree@all
    spns <- .getSltSpns(tree@ndlst)
    res <- .getNdsPrdstsFrmMtrx(
      tree@ndmtrx, all_ids, ids, spns,
      parallel, progress
    )
  } else {
    res <- .getNdsPrdstsFrmLst(tree@ndlst,
      prinds = tree@prinds,
      ids = ids, parallel, progress
    )
  }
  res
}

#' @name getNdsSlt
#' @title Get a node slot for multiple nodes
#' @description Returns the values of named slot as a vector for atomic values, else list.
#' @details Returned object depends on name, either character, vector or numeric. Parallelizable.
#' Default node slots are: id, spn, prid, ptid and txnym.
#' @param tree \code{TreeMan} object
#' @param slt_nm slot name
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdSlt}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsSlt(tree, slt_nm = "spn", ids = tree["tips"]) # return spans of all tips
getNdsSlt <- function(tree, slt_nm, ids, parallel = FALSE, progress = "none") {
  .get <- function(i) {
    getNdSlt(tree, slt_nm, ids[i])
  }
  l_data <- data.frame(i = 1:length(ids), stringsAsFactors = FALSE)
  res <- plyr::mlply(
    .data = l_data, .fun = .get, .parallel = parallel,
    .progress = progress
  )
  res <- res[1:length(res)]
  if (all(vapply(res, length, integer(1)) == 1)) {
    res <- unlist(res, recursive = FALSE)
  }
  names(res) <- ids
  res
}

#' @name getNdsKids
#' @title Get children IDs for multiple nodes
#' @description Return the node ids of all tips that descend from each node in \code{ids}.
#' @details Returns a list, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdKids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsKids(tree, id = tree["nds"])
getNdsKids <- function(tree, ids, parallel = FALSE,
                       progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    res <- .getNdsKidsFrmMtrx(
      tree@ndmtrx, tree@all,
      ids, tree@tips, parallel, progress
    )
  } else {
    res <- .getNdsKidsFrmLst(tree@ndlst,
      ids = ids,
      prinds = tree@prinds, tinds = tree@tinds,
      parallel = parallel, progress = progress
    )
  }
  res
}

#' @name getNdsAge
#' @title Get ages for multiple nodes
#' @description Return the age for \code{ids}.
#' @details Returns a vector, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdAge}},
#' \code{\link{getSpnAge}},
#' \code{\link{getSpnsAge}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsAge(tree, ids = tree["nds"], tree_age = getAge(tree))
getNdsAge <- function(tree, ids, tree_age,
                      parallel = FALSE,
                      progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    spns <- .getSltSpns(tree@ndlst)
    res <- .getNdsPrdstsFrmMtrx(
      tree@ndmtrx, tree@all,
      ids, spns, parallel, progress
    )
    res <- tree_age - res
  } else {
    res <- .getNdsPrdstsFrmLst(tree@ndlst,
      ids = ids,
      prinds = tree@prinds,
      parallel = parallel,
      progress = progress
    )
    res <- tree_age - res
  }
  res
}

#' @name getSpnsAge
#' @title Get age ranges for multiple nodes
#' @description Return start and end ages for \code{ids} from
#' when they first appear to when they split
#' @details Returns a dataframe, parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param tree_age numeric value of known age of tree
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdAge}},
#' \code{\link{getNdsAge}},
#' \code{\link{getSpnAge}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # all nodes but root
#' ids <- tree["nds"][tree["nds"] != tree["root"]]
#' getSpnsAge(tree, ids = ids, tree_age = getAge(tree))
getSpnsAge <- function(tree, ids, tree_age,
                       parallel = FALSE, progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    spns <- .getSltSpns(tree@ndlst)
    end <- .getNdsPrdstsFrmMtrx(
      tree@ndmtrx, tree@all,
      ids, spns, parallel, progress
    )
  } else {
    end <- .getNdsPrdstsFrmLst(tree@ndlst,
      ids = ids,
      prinds = tree@prinds,
      parallel = parallel,
      progress = progress
    )
  }
  spns <- getNdsSlt(tree, slt_nm = "spn", ids = ids, parallel)
  start <- end - spns
  end <- tree_age - end
  start <- tree_age - start
  data.frame(spn = ids, start, end, row.names = NULL)
}

#' @name getNdsPrids
#' @title Get pre-nodes for multiple nodes
#' @description Return node ids for connecting \code{id} to root.
#' @details Returns a list, parallizable. The function will work faster
#' if \code{ordrd} is FALSE.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param ordrd logical, ensure returned prids are ordered ID to root
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPrids}},
#' \code{\link{getNdPtids}},
#' \code{\link{getNdsPtids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' getNdsPrids(tree, ids = tree["tips"])
getNdsPrids <- function(tree, ids, ordrd = FALSE,
                        parallel = FALSE, progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1 & !ordrd) {
    res <- .getNdsPridsFrmMtrx(
      tree@ndmtrx, tree@all,
      ids, parallel, progress
    )
  } else {
    res <- .getNdsPridsFrmLst(tree@ndlst,
      ids = ids,
      prinds = tree@prinds, parallel = parallel,
      progress = progress
    )
  }
  res
}

#' @name getNdsPtids
#' @title Get post-nodes to tips for multiple nodes
#' @description Return node ids for connecting \code{ids} to kids.
#' @details Returns a list, parallizable.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{getNdPtids}},
#' \code{\link{getNdPrids}},
#' \code{\link{getNdsPrids}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # get all nodes to tip for all nodes
#' getNdsPtids(tree, ids = tree["nds"])
getNdsPtids <- function(tree, ids, parallel = FALSE, progress = "none") {
  if (!is.null(tree@ndmtrx) & length(ids) > 1) {
    res <- .getNdsPtidsFrmMtrx(
      tree@ndmtrx, tree@all,
      ids, parallel, progress
    )
  } else {
    res <- .getNdsPtidsFrmLst(tree@ndlst,
      ids = ids,
      prinds = tree@prinds, parallel = parallel,
      progress = progress
    )
  }
  res
}

# bipartitions
#' @name getBiprts
#' @title Get the sets of labels for each bipartition in tree
#' @description Returns a list of tip IDs for each branch in the tree. Options
#' allow the user to act as if the root is not present and to use a universal
#' code for comparing between trees.
#' @details Setting \code{root} to FALSE will ignore the bipartitions created by
#' the root. Setting \code{universal} to TRUE will return a vector of 0s and 1s,
#' not a list of tips. These codes will always begin with 1, and will allow for
#' the comparison of splits between trees as they do not have "chiralty", so to
#' speak.
#' @param tree \code{TreeMan} object
#' @param tips vector of tips IDs to use for bipartitions
#' @param root Include the root for the bipartitions? Default TRUE.
#' @param universal Create a code for comparing between trees
#' @seealso
#' \code{\link{calcDstRF}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # get all of the tip IDs for each branch in the rooted tree
#' (getBiprts(tree))
#' # ignore the root and get bipartitions for unrooted tree
#' (getBiprts(tree, root = FALSE))
#' # use the universal code for comparing splits between trees
#' (getBiprts(tree, root = FALSE, universal = TRUE))
getBiprts <- function(tree, tips = tree@tips, root = TRUE, universal = FALSE) {
  kids <- getNdsKids(tree = tree, ids = tree@nds)
  res <- lapply(X = kids, FUN = function(x) tips %in% x)
  if (!root) {
    n <- vapply(X = res, FUN = sum, FUN.VALUE = integer(1))
    # drop splits consisting of all tips or just 1
    res <- res[n < (length(tips) - 1) & n > 1]
  }
  if (universal) {
    res <- unname(vapply(X = res, FUN = function(x) {
      if (!x[[1]]) {
        x <- !x
      }
      paste0(as.integer(x), collapse = "")
    }, FUN.VALUE = character(1)))
    res <- unique(res)
  } else {
    res <- lapply(X = res, FUN = function(x) {
      list(tree@tips[x], tree@tips[!x])
    })
  }
  res
}

# ULTRAMETRIC
#' @name isUltrmtrc
#' @title Is tree ultrametric?
#' @description Return TRUE if all tips end at 0, else FALSE.
#' @details Returns a boolean. This function works in the background
#' for the \code{['ultr']} slot in a \code{TreeMan} object.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getLvng}}, \code{\link{getDcsd}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' (isUltrmtrc(tree))
isUltrmtrc <- function(tree, tol = 1e-8) {
  dead <- .livingOrDesceased(tree, tol = tol, bool = FALSE)
  if (length(dead) > 0) {
    return(FALSE)
  }
  TRUE
}

# EXTINCT/EXTANT
.livingOrDesceased <- function(tree, tol = 1e-8, bool) {
  if (!is.null(tree@ndmtrx)) {
    spns <- getNdsSlt(tree, "spn", names(tree@ndlst))
    tip_prdsts <- .getNdsPrdstsFrmMtrx(tree@ndmtrx, tree@all,
      tree@tips, spns,
      parallel = FALSE,
      progress = "none"
    )
  } else {
    tip_prdsts <- .getNdsPrdstsFrmLst(tree@ndlst, tree@tips,
      tree@prinds,
      parallel = FALSE,
      progress = "none"
    )
  }
  age <- max(tip_prdsts)
  extant_is <- (age - tip_prdsts) <= tol
  living <- names(extant_is)[extant_is]
  deceased <- tree@tips[!tree@tips %in% living]
  if (bool) {
    return(living)
  }
  deceased
}

#' @name getDcsd
#' @title Get extinct tips from a tree
#' @description Return all extinct tip \code{ID}s.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getLvng}}, \code{\link{isUltrmtrc}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' (getDcsd(tree))
getDcsd <- function(tree, tol = 1e-8) {
  .livingOrDesceased(tree = tree, tol = tol, bool = FALSE)
}

#' @name getLvng
#' @title Get extant tips from a tree
#' @description Return all extant tip \code{ID}s.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tol zero tolerance
#' @seealso
#' \code{\link{getDcsd}}, \code{\link{isUltrmtrc}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' (getLvng(tree))
getLvng <- function(tree, tol = 1e-8) {
  .livingOrDesceased(tree = tree, tol = tol, bool = TRUE)
}

# SINGLE ND
# TODO: bring outgroup, parent and path into terminological line with getNd(s)

#' @name getPrnt
#' @title Get parent
#' @description Return parental (most recent common ancestor) node id for \code{ids}.
#' @details Returns a character.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \code{\link{getSubtree}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # choosing ids from the two main branches of apes allows to find the parent for all apes
#' ape_id <- getPrnt(mammals, ids = c("Homo_sapiens", "Hylobates_concolor"))
getPrnt <- function(tree, ids) {
  # using ndlst guarrantees order
  prids <- .getNdsPridsFrmLst(tree@ndlst,
    ids = ids, prinds = tree@prinds,
    parallel = FALSE, progress = "none"
  )
  rf <- prids[[1]]
  mn_rnk <- 0
  for (n in prids[-1]) {
    rnk <- min(match(n, rf), na.rm = TRUE)
    if (rnk > mn_rnk) mn_rnk <- rnk
  }
  rf[mn_rnk]
}

#' @name getPath
#' @title Get path between nodes
#' @description Return node ids for connecting \code{from} to \code{to}.
#' @details Returns a vector, first id is \code{from} to \code{to}.
#' @param tree \code{TreeMan} object
#' @param from starting node id
#' @param to ending node id
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # what's the phylogenetic distance from humans to gorillas?
#' ape_id <- getPrnt(mammals, ids = c("Homo_sapiens", "Hylobates_concolor"))
#' pth <- getPath(mammals, from = "Homo_sapiens", to = "Gorilla_gorilla")
#' sum(getNdsSlt(mammals, ids = pth, slt_nm = "spn"))
getPath <- function(tree, from, to) {
  pre_1 <- c(from, getNdPrids(tree, from))
  pre_2 <- c(to, getNdPrids(tree, to))
  parent <- pre_1[which(pre_1 %in% pre_2)[1]]
  path_1 <- pre_1[!pre_1 %in% pre_2]
  path_2 <- pre_2[!pre_2 %in% pre_1]
  path_2 <- path_2[length(path_2):1]
  c(path_1, parent, path_2)
}

#' @name getOtgrp
#' @title Get outgroup
#' @description Return the outgroup based on a tree and a vector of IDs.
#' @details Returns a id, character. If there are multiple possible outgroups, returns NULL.
#' @param tree \code{TreeMan} object
#' @param ids vector of node ids
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # orangutan is an outgroup wrt humans and chimps
#' getOtgrp(mammals, ids = c("Homo_sapiens", "Pan_troglodytes", "Pongo_pygmaeus"))
getOtgrp <- function(tree, ids) {
  .cntr <- function(id) {
    kids <- getNdKids(tree, id)
    sum(ids %in% kids)
  }
  prnt <- getPrnt(tree, ids)
  ptids <- tree@ndlst[[prnt]][["ptid"]]
  cnts <- vapply(ptids, .cntr, integer(1))
  outnd <- names(cnts)[which.min(cnts)]
  kids <- getNdKids(tree, outnd)
  if (length(kids) == 0) {
    return(outnd)
  }
  outgroup <- ids[ids %in% kids]
  if (length(outgroup) > 1) {
    return(NULL)
  }
  outgroup
}

# SPECIAL

#' @name getUnqNds
#' @title Get unique nodes represented by tips
#' @description Return a list of IDs for any node that are represented by tip IDs given.
#' @details Returns a vector.
#' @param tree \code{TreeMan} object
#' @param tids vector of tip IDs
#' @seealso
#' \code{\link{getCnnctdNds}}, \code{\link{calcFrPrp}},
#' \code{\link{calcPhyDv}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' unqnds <- getUnqNds(tree, c("t1", "t2"))
getUnqNds <- function(tree, tids) {
  rmng <- tree@tips[!tree@tips %in% tids]
  ignr <- c(unique(unlist(getNdsPrids(tree, tids))), tids)
  rmng <- c(unique(unlist(getNdsPrids(tree, rmng))), rmng)
  ignr[!ignr %in% rmng]
}

#' @name getCnnctdNds
#' @title Get all nodes connected by given tips
#' @description Return a vector of IDs of all nodes that are connected to tip IDs given.
#' @details Returns a vector. This function is the basis for \code{calcPhyDv()}, it determines
#' the unique set of nodes connected for a set of tips.
#' @param tree \code{TreeMan} object
#' @param tids vector of tip IDs
#' @seealso
#' \code{\link{getUnqNds}}, \code{\link{calcFrPrp}},
#' \code{\link{calcPhyDv}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' cnntdnds <- getCnnctdNds(tree, c("t1", "t2"))
getCnnctdNds <- function(tree, tids) {
  prids <- c(
    unlist(getNdsPrids(tree, tids)),
    tids
  )
  counts <- table(prids)
  names(counts)[counts < length(tids)]
}

#' @name getNdsFrmTxnyms
#' @title Get IDs for nodes represented txnyms
#' @description Return a list of IDs for any node that contains the given txnyms.
#' @details Returns a list. Txnyms must be spelt correctly.
#' @param tree \code{TreeMan} object
#' @param txnyms vector of taxonomic group names
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{setTxnyms}}, \code{\link{searchTxnyms}},
#' \code{\link{getNdsLng}}, \code{\link{getNdLng}}
#' @export
#' @examples
#'
#' data(mammals)
#' # what ID represents the apes?
#' getNdsFrmTxnyms(mammals, "Hominoidea")
getNdsFrmTxnyms <- function(tree, txnyms) {
  # get nd id(s) for taxonyms
  .get <- function(id, txnym, ...) {
    for (t in txnyms) {
      if (t %in% txnym) {
        res[[t]] <<- c(res[[t]], id)
      }
    }
  }
  res <- list()
  plyr::m_ply(tree@ndlst, .fun = .get)
  res
}

#' @name getSubtree
#' @title Get subtree
#' @description Return tree descending from \code{id}.
#' @details Returns a \code{TreeMan}, parallelizable. \code{id} must be an internal node.
#' @param tree \code{TreeMan} object
#' @param id node id
#' @seealso
#' \code{\link{getPrnt}}, \code{\link{addClade}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # get tree of apes
#' ape_id <- getPrnt(mammals, ids = c("Homo_sapiens", "Hylobates_concolor"))
#' apes <- getSubtree(mammals, id = ape_id)
#' summary(apes)
getSubtree <- function(tree, id) {
  if (!id %in% tree@nds) {
    stop("`id` is not an internal node")
  }
  ids <- c(id, getNdPtids(tree, id))
  ndlst <- tree@ndlst[ids]
  ndlst[[id]][["prid"]] <- id
  ndlst[[id]][["spn"]] <- 0
  new_tree <- new("TreeMan",
    ndlst = ndlst, root = id,
    ndmtrx = NULL
  )
  new_tree <- pstMnp(new_tree)
  new_tree <- updateSlts(new_tree)
  new_tree
}

# TREE FUNCTIONS

#' @name getAge
#' @title Get age of tree
#' @description Returns age, numeric, of tree
#' @details Calculates the age of a tree, determined as the maximum tip to root
#' distance.
#' @param tree \code{TreeMan} object
#' @param parallel logical, make parallel?
#' @seealso
#' \code{\link{updateSlts}},
#' \url{https://github.com/DomBennett/treeman/wiki/get-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' (getAge(tree))
getAge <- function(tree, parallel = FALSE) {
  tids <- tree@tips
  if (!is.null(tree@ndmtrx)) {
    all_ids <- tree@all
    spns <- .getSltSpns(tree@ndlst)
    res <- .getTreeAgeFrmMtrx(tree@ndmtrx, all_ids, tids, spns, parallel)
  } else {
    res <- .getTreeAgeFrmLst(tree@ndlst,
      prinds = tree@prinds,
      tids = tids, parallel
    )
  }
  res
}


#' @name ultrTree
#' @title Make tree ultrametric
#' @description Returns a tree with all tips ending at time 0
#' @details Re-calculates the branch lengths in the tree so that all
#' tips are brought to the same time point: all species are extant.
#' @param tree \code{TreeMan} object
#' @seealso
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' (getDcsd(tree)) # list all extinct tips
#' tree <- ultrTree(tree)
#' (getDcsd(tree)) # list all extinct tips
ultrTree <- function(tree) {
  # bring tips to maximum possible length
  tip_dpths <- vapply(getNdsPrids(tree, tree@tips), length, integer(1))
  tip_spns <- (tree["ntips"] - tip_dpths)
  nd_spns <- rep(1, tree@nnds)
  names(nd_spns) <- tree@nds
  spns <- c(tip_spns, nd_spns)
  tree <- setNdsSpn(tree, ids = names(spns), vals = spns)
  updateSlts(tree)
}

#' @name rmNodes
#' @title Remove nodes from a tree
#' @description Returns a tree with a node ID(s) removed
#' @details Removes nodes in a tree. Joins the nodes following to
#' the nodes preceding the node to be removed. Creates polytomies.
#' Warning: do not use this function to remove tip nodes, this create a
#' corrupted tree.
#' @param tree \code{TreeMan} object
#' @param nids internal node IDs
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmTips}},
#'  \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- rmNodes(tree, "n3")
#' summary(tree) # tree is now polytmous
rmNodes <- function(tree, nids, progress = "none") {
  .rmNode <- function(nid) {
    ptids <- ndlst[[nid]][["ptid"]]
    prid <- ndlst[[nid]][["prid"]]
    if (tree@wspn) {
      for (ptid in ptids) {
        ndlst[[ptid]][["spn"]] <-
          ndlst[[ptid]][["spn"]] +
          ndlst[[nid]][["spn"]]
      }
    }
    for (ptid in ptids) {
      ndlst[[ptid]][["prid"]] <- prid
    }
    new_ptids <- ndlst[[prid]][["ptid"]]
    new_ptids <- new_ptids[new_ptids != nid]
    new_ptids <- c(new_ptids, ptids)
    ndlst[[prid]][["ptid"]] <- new_ptids
    ndlst <<- ndlst[names(ndlst) != nid]
  }
  if (tree@root %in% nids) {
    stop("Cannot remove root.")
  }
  ndlst <- tree@ndlst
  plyr::m_ply(.data = nids, .fun = .rmNode, .progress = progress)
  bool <- tree@all %in% names(ndlst)
  tree@ndlst <- ndlst
  tree <- pstMnp(tree)
  tree <- updateSlts(tree)
  if (!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  tree
}

#' @name rmTips
#' @title Remove tips from a tree
#' @description Returns a tree with a tip ID(s) removed
#' @details Removes tips in a tree. Set drp_intrnl to FALSE to convert
#' internal nodes into new tips. Warning: do not use this function to remove
#' internal nodes, this create a corrupted tree.
#' @param tree \code{TreeMan} object
#' @param tids tip IDs
#' @param drp_intrnl Boolean, drop internal branches, default FALSE
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmNodes}},
#'  \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- rmTips(tree, "t1")
#' summary(tree)
#' # running the function using an internal
#' # node will create a corrupted tree
#' tree <- rmTips(tree, "n3")
#' # run summary() to make sure a change has
#' # not created a corruption
#' # summary(tree)
rmTips <- function(tree, tids, drp_intrnl = TRUE, progress = "none") {
  # internal
  .rmTip <- function(tid) {
    # get sister IDs
    sids <- .getNdSstrFrmLst(ndlst, tid)
    # get prid
    prid <- ndlst[[tid]][["prid"]][[1]]
    # remove tid
    ndlst <<- ndlst[names(ndlst) != tid]
    ndlst[[prid]][["ptid"]] <<-
      ndlst[[prid]][["ptid"]][ndlst[[prid]][["ptid"]] != tid]
    # remove prnd if specified and not polytomous
    if (drp_intrnl & length(sids) == 1) {
      ptid <- ndlst[[prid]][["ptid"]][[1]]
      ndlst[[ptid]][["spn"]] <<- ndlst[[prid]][["spn"]] +
        ndlst[[ptid]][["spn"]]
      if (prid != rid) {
        gprid <- ndlst[[prid]][["prid"]][[1]]
        ndlst[[ptid]][["prid"]] <<- gprid
        g_ptids <- ndlst[[gprid]][["ptid"]]
        g_ptids <- g_ptids[g_ptids != prid]
        ndlst[[gprid]][["ptid"]] <<- c(g_ptids, ptid)
      } else {
        # if prid to be dropped is root, set root to ptid
        rid <<- ptid
        ndlst[[ptid]][["prid"]] <<- ptid
      }
      ndlst <<- ndlst[names(ndlst) != prid]
    }
  }
  ndlst <- tree@ndlst
  rid <- tree@root
  plyr::m_ply(.data = tids, .fun = .rmTip, .progress = progress)
  bool <- tree@all %in% names(ndlst)
  tree@ndlst <- ndlst
  tree@root <- rid
  tree <- pstMnp(tree)
  tree <- updateSlts(tree)
  if (!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  tree
}

#' @name addTip
#' @title Add tip to a tree
#' @description Returns a tree with a new tip ID added
#' @details User must provide new tip ID, the ID of the node
#' which will become the new tip's sister, and new branch lengths.
#' The tip ID must only contain letters numbers and underscores.
#' Optionally, user can specify the IDs for the new parental internal nodes.
#' Ensure that the \code{strt_age} is greater than the \code{end_age}, and that
#' the \code{strt_age} falls within the age span of the sister ID. Otherwise, negative
#' spns may be produced leading to an error.
#' Note, returned tree will not have a node matrix.
#' Note, providing negative end ages will increase the age of the tree.
#' @param tree \code{TreeMan} object
#' @param tid tip ID
#' @param sid ID of node that will become new tip sisters
#' @param strt_age timepoint at which new tips first appear in the tree
#' @param end_age timepoint at which new tips end appear in the tree, default 0.
#' @param tree_age age of tree
#' @param pid parent ID (default is 'p_' + tid)
#' @seealso
#' \code{\link{rmTips}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree_age <- getAge(tree)
#' possible_ages <- getSpnAge(tree, "t1", tree_age)
#' start_age <- runif(1, possible_ages[["end"]], possible_ages[["start"]])
#' end_age <- possible_ages[["end"]]
#' tree <- addTip(tree,
#'   tid = "t11", sid = "t1", strt_age = start_age,
#'   end_age = end_age, tree_age = tree_age
#' )
#' summary(tree)
addTip <- function(tree, tid, sid, strt_age = NULL,
                   end_age = 0, tree_age = NULL,
                   pid = paste0("p_", tid)) {
  if (!is.numeric(strt_age) & tree@wspn) {
    stop("Valid strt_age not provided")
  }
  ndlst <- tree@ndlst
  # terminology
  # snd, sid -- old sister node and id
  # tnd, tid -- new tip node and id
  # pnd, pid -- new parent node and id
  # gpnd, gpid -- grand parent (prid of old sister)
  if (grepl("[^a-zA-Z_0-9]", tid)) {
    stop(paste0("Unsuitable characters in tid [", tid, "]"))
  }
  # init new nodes
  tnd <- list("id" = tid, "prid" = pid, "ptid" = character(), "spn" = 0)
  snd <- ndlst[[sid]]
  gpid <- snd[["prid"]][[1]]
  gpnd <- ndlst[[gpid]]
  pnd <- list("id" = pid, "spn" = 0)
  # update ptid
  gpnd[["ptid"]] <- gpnd[["ptid"]][!gpnd[["ptid"]] %in% snd[["id"]]]
  gpnd[["ptid"]] <- c(gpnd[["ptid"]], pnd[["id"]])
  pnd[["ptid"]] <- c(tid, sid)
  # set prid
  pnd[["prid"]] <- snd[["prid"]]
  snd[["prid"]] <- pid
  # add to ndlst
  ndlst[[tid]] <- tnd
  ndlst[[pid]] <- pnd
  ndlst[[sid]] <- snd
  ndlst[[gpid]] <- gpnd
  if (tree@wspn) {
    sspn <- ndlst[[sid]][["spn"]]
    prid <- ndlst[[tid]][["prid"]]
    gprid <- ndlst[[prid]][["prid"]]
    # calc
    gp_age <- getNdAge(tree, gprid, tree_age)
    tspn <- strt_age - end_age
    pspn <- gp_age - strt_age
    new_spsn <- abs(sspn - pspn)
    if (tspn < 0 | pspn < 0 | new_spsn < 0) {
      stop("Invalid ages given: negative spns")
    }
    ndlst[[tid]][["spn"]] <- tspn
    ndlst[[prid]][["spn"]] <- pspn
    ndlst[[sid]][["spn"]] <- new_spsn
  }
  tree@ndlst <- ndlst
  tree <- pstMnp(tree)
  tree <- rmNdmtrx(tree)
  tree <- updateSlts(tree)
  tree
}

#' @name addClade
#' @title Add clade to tree
#' @description Returns a tree with added clade
#' @details Add a \code{TreeMan} object to an existing \code{TreeMan}
#' object by specifying an ID at which to attach. If the id specified
#' is an internal node, then the original clade descending from that
#' node will be replaced. Before running, ensure no IDs are shared
#' between the \code{tree} and the \code{clade}, except for the IDs in the clade
#' of that tree that will be replaced.
#' Note, returned tree will not have a node matrix.
#' @param tree \code{TreeMan} object
#' @param id tip/node ID in tree to which the clade will be added
#' @param clade \code{TreeMan} object
#' @seealso
#' \code{\link{rmClade}}, \code{\link{getSubtree}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' t1 <- randTree(100)
#' # extract a clade
#' cld <- getSubtree(t1, "n2")
#' # remove the same clade
#' t2 <- rmClade(t1, "n2")
#' # add the clade again
#' t3 <- addClade(t2, "n2", cld)
#' # t1 and t3 should be the same
#' # note there is no need to remove a clade before adding
#' t3 <- addClade(t1, "n2", cld) # same tree
addClade <- function(tree, id, clade) {
  if (!id %in% tree@tips) {
    tree <- rmClade(tree, id)
  }
  cld_ids <- names(clade@ndlst)[names(clade@ndlst) != clade@root]
  if (any(names(tree@ndlst) %in% cld_ids) &
    any(cld_ids %in% names(tree@ndlst))) {
    stop(
      "IDs in `clade` exist in parts of `tree` not to be replaced.",
      " Consider checking for duplicates or renaming IDs in either `tree` or `clade`"
    )
  }
  cld_ptids <- clade@ndlst[[clade@root]][["ptid"]]
  for (cld_ptid in cld_ptids) {
    clade@ndlst[[cld_ptid]][["prid"]] <- id
  }
  cld_ndlst <- clade@ndlst[names(clade@ndlst) != clade@root]
  tree@ndlst[[id]][["ptid"]] <- cld_ptids
  tree@ndlst <- c(tree@ndlst, cld_ndlst)
  tree <- pstMnp(tree)
  tree <- rmNdmtrx(tree)
  updateSlts(tree)
}

#' @name rmClade
#' @title Remove a clade from a tree
#' @description Returns a tree with a clade removed
#' @details Inverse function of \code{getSubtree()}. Takes a tree
#' and removes a clade based on an internal node specified. Node
#' is specified with \code{id}, all descending nodes and tips are removed.
#' The resulting tree will replace the missing clade with a tip of \code{id}.
#' @param tree \code{TreeMan} object
#' @param id node ID parent of clade to be removed
#' @seealso
#' \code{\link{addClade}}, \code{\link{getSubtree}}, \code{\link{rmTips}}
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#'
#' t1 <- randTree(100)
#' # remove a clade
#' t2 <- rmClade(t1, "n2")
#' summary(t1)
#' summary(t2)
rmClade <- function(tree, id) {
  ptids <- getNdPtids(tree, id)
  bool <- !tree@all %in% ptids
  tree@ndlst <- tree@ndlst[!names(tree@ndlst) %in% ptids]
  tree@ndlst[[id]][["ptid"]] <- character()
  tree <- pstMnp(tree)
  if (!is.null(tree@ndmtrx)) {
    tree@ndmtrx <- bigmemory::as.big.matrix(tree@ndmtrx[bool, bool])
  }
  updateSlts(tree)
}

#' @name pinTips
#' @title Pin tips to a tree
#' @description Returns a tree with new tips added based on given lineages and time points
#' @details User must provide a vector of new tip IDs, a list of the ranked lineages
#' for these IDs (in ascending order) and a vector of end time points for each new ID
#' (0s for extant tips). The function expects the given tree to be taxonomically informed;
#' the \code{txnym} slot for every node should have a taxonomic label. The function takes
#' the lineage and tries to randomly add the new tip at the lowest point in the taxonomic rank
#' before the end time point. Note, returned tree will not have a node matrix.
#' @param tree \code{TreeMan} object
#' @param tids new tip ids
#' @param lngs list of vectors of the lineages of each tid (ordered high to low rank)
#' @param end_ages end time points for each tid
#' @param tree_age age of tree
#' @seealso
#' \code{\link{addTip}}, \code{\link{rmTips}},
#' \url{https://github.com/DomBennett/treeman/wiki/manip-methods}
#' @export
#' @examples
#' # see https://github.com/DomBennett/treeman/wiki/Pinning-tips for a detailed example
pinTips <- function(tree, tids, lngs, end_ages, tree_age) {
  .getPtntls <- function(lng, end) {
    sccs <- FALSE
    # loop backwards through taxonomy
    # genus --> family --> ....
    for (i in length(lng):1) {
      pull <- txnyms %in% lng[i]
      if (sum(pull) == 0) {
        next
      }
      ptntls <- names(txnyms)[pull]
      if (length(ptntls) == 1) {
        prnt <- ptntls
      } else {
        # assumes monophylly
        prnt <- ptntls[which.max(spn_data[ptntls, "end"])]
      }
      ptntls <- c(prnt, getNdPtids(tree, prnt))
      ptntls <- ptntls[ptntls != rid]
      pull <- spn_data[ptntls, "start"] > end
      if (any(pull)) {
        ptntls <- ptntls[pull]
        sccs <- TRUE
        break
      }
    }
    if (sccs) {
      return(ptntls)
    }
    NULL
  }
  .getPTxnym <- function(tip_txnym, sid) {
    gp_txnym <- txnyms[[getNdSlt(tree, "prid", sid)]]
    s_txnym <- txnyms[[sid]]
    if (s_txnym == tip_txnym) {
      pid_txnym <- tip_txnym
    } else {
      pid_txnym <- gp_txnym
    }
    pid_txnym
  }
  .pin <- function(i) {
    tid <- tids[i]
    end <- end_ages[i]
    lng <- lngs[[i]]
    ptntls <- .getPtntls(lng, end)
    if (is.null(ptntls)) {
      message(paste0("[", tid, "] could not be added"))
      return(NULL)
    }
    rngs <- spn_data[ptntls, , drop = FALSE]
    rngs[rngs[, "end"] <= end, "end"] <- end
    # pinning is based on branch length
    # this is not a model, it just ensures
    # taxonomically matching branch lengths
    # of the tree have equal chance.
    prbs <- rngs[, "start"] - rngs[, "end"]
    if (sum(prbs) < 10 - 8) {
      # cannot add where no branch is available, prbs eqv. to branch length
      message(paste0("[", tid, "] could not be added"))
      return(NULL)
    }
    sid <- as.vector(sample(ptntls, prob = prbs, size = 1))
    start <- runif(min = rngs[sid, "end"], max = rngs[sid, "start"], n = 1)
    # taxnomy of tip and parent tip based grandparent
    tip_txnym <- lng[length(lng)]
    pid_txnym <- .getPTxnym(tip_txnym, sid)
    pid <- paste0("p_", tid, sep = "")
    # add tip
    tree <- addTip(tree,
      tid = tid, sid = sid, strt_age = start,
      end_age = end, pid = pid, tree_age = tree_age
    )
    tree@ndlst[[tid]][["txnym"]] <- tip_txnym
    tree@ndlst[[pid]][["txnym"]] <- pid_txnym
    # add to spn_data
    tid_spn <- getSpnAge(tree, tid, tree_age)
    spn_data[tid, "start"] <<- tid_spn[, "start"]
    spn_data[tid, "end"] <<- tid_spn[, "end"]
    pid_spn <- getSpnAge(tree, pid, tree_age)
    spn_data[pid, "start"] <<- pid_spn[, "start"]
    spn_data[pid, "end"] <<- pid_spn[, "end"]
    sid_spn <- getSpnAge(tree, sid, tree_age)
    spn_data[sid, "start"] <<- sid_spn[, "start"]
    spn_data[sid, "end"] <<- sid_spn[, "end"]
    # add to txnyms list
    txnyms[[tid]] <<- tip_txnym
    txnyms[[pid]] <<- pid_txnym
    # push tree out
    tree <<- tree
  }
  .testLngs <- function(lng) {
    for (l in lng) {
      if (grepl("[^a-zA-Z_0-9]", l)) {
        stop(paste0("Unsuitable characters in [", l, "]"))
      }
    }
    NULL
  }
  if (!tree@wtxnyms) {
    stop("tree has no txnyms")
  }
  if (any(end_ages < 0)) {
    warning("One or more end ages are less than zero, this will change the age of tree.")
  }
  # make sure lineages are right
  mapply(.testLngs, lngs)
  # generate taxonomy and span data
  txnyms <- getNdsSlt(tree, "txnym", tree@all)
  txnyms <- c(txnyms, rep(NA, length(tids) * 2))
  names(txnyms) <- c(names(tree@ndlst), tids, paste0("p_", tids))
  spn_data <- matrix(NA,
    nrow = (length(tids) * 2) + tree@nall,
    ncol = 2
  )
  colnames(spn_data) <- c("start", "end")
  tmp_spn_data <- getSpnsAge(tree, tree@all, tree_age)
  rownames(spn_data) <- c(tree@all, tids, paste0("p_", tids))
  spn_data[tree@all, "start"] <- tmp_spn_data[["start"]]
  spn_data[tree@all, "end"] <- tmp_spn_data[["end"]]
  rm(tmp_spn_data)
  rid <- tree@root
  # add oldest to youngest
  ordrd <- order(end_ages, decreasing = TRUE)
  plyr::m_ply(ordrd, .pin)
  tree
}

# prevent sudden crashes by stopping incorrect arguments going to C code
.preCtest <- function(...) {
  argg <- c(as.list(environment()), list(...))
  bool <- FALSE
  for (arg in argg) {
    if (!is.numeric(arg) | length(arg) < 1) {
      bool <- TRUE
    }
  }
  if (bool) {
    stop("1 or more arguments are of the wrong type")
  }
  NULL
}

# MULTIPLE NDS

.getNdsPtidsFrmLst <- function(ndlst, ids, prinds, parallel, progress) {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  out <- plyr::mlply(
    .data = l_data, .fun = .getNdPtidsFrmLst, ndlst = ndlst,
    prinds = prinds, .parallel = parallel, .progress = progress
  )
  names(out) <- attr(out, "split_labels")[, 1]
  res <- out[1:length(out)]
  res
}

.getNdsPridsFrmLst <- function(ndlst, ids, prinds,
                               parallel, progress) {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  out <- plyr::mlply(
    .data = l_data, .fun = .getNdPridsFrmLst, ndlst = ndlst,
    prinds = prinds, .parallel = parallel, .progress = progress
  )
  names(out) <- attr(out, "split_labels")[, 1]
  res <- out[1:length(out)]
  res
}

.getNdsPDFrmLst <- function(ndlst, ids, prinds,
                            parallel, progress) {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  out <- plyr::mdply(
    .data = l_data, .fun = .getNdPDFrmLst, prinds = prinds,
    ndlst = ndlst, .parallel = parallel, .progress = progress
  )
  res <- out[, 2]
  names(res) <- out[, 1]
  res
}

.getNdsKidsFrmLst <- function(ndlst, ids, prinds, tinds,
                              parallel, progress) {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  res <- plyr::mlply(
    .data = l_data, .fun = .getNdKidsFrmLst,
    ndlst = ndlst, prinds = prinds,
    tinds = tinds, .parallel = parallel, .progress = progress
  )
  names(res) <- ids
  res[1:length(res)]
}

.getNdsPrdstsFrmLst <- function(ndlst, ids, prinds,
                                parallel, progress) {
  l_data <- data.frame(id = ids, stringsAsFactors = FALSE)
  out <- plyr::mdply(
    .data = l_data, .fun = .getNdPrdstsFrmLst, prinds = prinds,
    ndlst = ndlst, .parallel = parallel, .progress = progress
  )
  res <- out[, 2]
  names(res) <- out[, 1]
  res
}

# SINGLE ND

.getNdSstrFrmLst <- function(ndlst, id) {
  prid <- ndlst[[id]][["prid"]][[1]]
  ptids <- ndlst[[prid]][["ptid"]]
  ptids[ptids != id]
}

#' @useDynLib phylotaR cGetNdPrids
.getNdPridsFrmLst <- function(ndlst, prinds, id) {
  prid <- ndlst[[id]][["prid"]]
  nids <- names(ndlst)
  prind <- which(nids == prid)
  .preCtest(prind, prinds)
  res <- .Call("cGetNdPrids",
    PACKAGE = "phylotaR",
    as.integer(prind),
    as.integer(prinds)
  )
  nids[res]
}

.getNdPrdstsFrmLst <- function(ndlst, prinds, id) {
  prids <- .getNdPridsFrmLst(ndlst, prinds, id)
  sum(vapply(ndlst[prids], function(x) x[["spn"]], numeric(1))) +
    ndlst[[id]][["spn"]]
}

#' @useDynLib phylotaR cGetNdPtids
.getNdPtidsFrmLst <- function(ndlst, prinds, id) {
  nids <- names(ndlst)
  id <- which(nids == id)
  .preCtest(id, prinds)
  res <- .Call("cGetNdPtids",
    PACKAGE = "phylotaR",
    as.integer(id),
    as.integer(prinds)
  )
  nids[which(res > 0)]
}

.getNdKidsFrmLst <- function(ndlst, prinds, tinds, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prinds, id)
  tids <- names(ndlst)[tinds]
  ptids[ptids %in% tids]
}

.getNdPDFrmLst <- function(ndlst, prinds, id) {
  ptids <- .getNdPtidsFrmLst(ndlst, prinds, id)
  if (length(ptids) > 0) {
    res <- sum(vapply(ndlst[ptids], function(x) x[["spn"]], numeric(1)))
  } else {
    res <- 0
  }
  res
}

# TREE FUNCTIONS

.getTreeAgeFrmLst <- function(ndlst, prinds, tids, parallel) {
  prdsts <- .getNdsPrdstsFrmLst(ndlst,
    ids = tids, prinds = prinds,
    parallel = parallel, progress = "none"
  )
  max(prdsts)
}

# SLOT

.getSltSpns <- function(ndlst) {
  .get <- function(x) {
    x[["spn"]]
  }
  res <- vapply(ndlst, .get, numeric(1))
  names(res) <- NULL
  res
}

# INDS

.getPrinds <- function(ndlst) {
  # pre-node index
  .get <- function(x) {
    x[["prid"]]
  }
  res <- vapply(ndlst, .get, character(1))
  match(res, names(ndlst))
}

.getTinds <- function(ndlst) {
  # tip-node index
  .get <- function(x) {
    length(x[["ptid"]]) == 0
  }
  res <- vapply(ndlst, .get, logical(1))
  names(res) <- NULL
  which(res)
}

# SPECIAL

#' @useDynLib phylotaR cGetNdmtrx
.getNdmtrxFrmLst <- function(ndlst, shared = FALSE, ...) {
  # return matrix of 01s for ids that descend from
  message(
    "Note, trees with `ndmtrx` cannot be saved and loaded using `save()` or `savehistory()`.",
    " Loading from these files may cause unusual behaviour."
  )
  prids <- vapply(ndlst, function(x) x[["prid"]], character(1))
  nids <- names(prids)
  prids <- match(prids, nids)
  qry_ids <- 1:length(nids)
  .preCtest(length(nids), qry_ids, prids)
  res <- .Call("cGetNdmtrx",
    PACKAGE = "phylotaR",
    as.integer(length(nids)),
    as.integer(qry_ids),
    as.integer(prids)
  )
  res <- bigmemory::as.big.matrix(res, shared = shared, ...)
  res
}
# Attemp for making getNdsMat run in parallel
# ... actually made it slower
# ntids <- length(tids)
# n <- foreach::getDoParWorkers()
# nparts <- ntids %/% n
# parts <- c(seq(1, ntids - 1, nparts), ntids + 1)
# res <- foreach (i=2:length(parts), .combine="cbind") %dopar% {
#   tids <- tids[parts[i-1]:(parts[i] - 1)]
#   res <- .Call("cGetNdsMat", PACKAGE="phylotaR",
#                as.integer(length(nids)),
#                as.integer(tids),
#                as.integer(prids))
#   res
# }


# MULTIPLE NDS
.getNdsPrdstsFrmMtrx <- function(ndmtrx, all_ids, ids, spns,
                                 parallel, progress) {
  .get <- function(x) {
    sum(spns[as.logical(x)])
  }
  res <- plyr::adply(ndmtrx[, all_ids %in% ids],
    .margins = 2,
    .fun = .get, .parallel = parallel, .progress = progress
  )[, 2]
  res <- res + spns[all_ids %in% ids]
  names(res) <- all_ids[all_ids %in% ids]
  res
}

.getNdsKidsFrmMtrx <- function(ndmtrx, all_ids, ids, tids,
                               parallel, progress) {
  .get <- function(x) {
    all_ids[as.logical(x) & all_ids %in% tids]
  }
  res <- plyr::alply(ndmtrx[all_ids %in% ids, ],
    .margins = 1,
    .fun = .get, .parallel = parallel, .progress = progress
  )
  names(res) <- all_ids[all_ids %in% ids]
  res <- res[1:length(res)]
  res
}

.getNdsPtidsFrmMtrx <- function(ndmtrx, all_ids, ids,
                                parallel, progress) {
  .get <- function(x) {
    all_ids[as.logical(x)]
  }
  res <- plyr::alply(ndmtrx[all_ids %in% ids, ],
    .margins = 1,
    .fun = .get, .parallel = parallel, .progress = progress
  )
  names(res) <- all_ids[all_ids %in% ids]
  res <- res[1:length(res)]
  res
}

.getNdsPridsFrmMtrx <- function(ndmtrx, all_ids, ids,
                                parallel, progress) {
  .get <- function(x) {
    all_ids[as.logical(x)]
  }
  res <- plyr::alply(ndmtrx[, all_ids %in% ids],
    .margins = 2,
    .fun = .get, .parallel = parallel, .progress = progress
  )
  names(res) <- all_ids[all_ids %in% ids]
  res <- res[1:length(res)]
  res
}

.getNdsPDFrmMtrx <- function(ndmtrx, all_ids, ids, spns,
                             parallel, progress) {
  .get <- function(x) {
    sum(spns[as.logical(x)])
  }
  res <- plyr::adply(ndmtrx[all_ids %in% ids, ],
    .margins = 1,
    .fun = .get, .parallel = parallel, .progress = progress
  )[, 2]
  names(res) <- all_ids[all_ids %in% ids]
  res <- res[1:length(res)]
  res
}

# TREE

.getTreeAgeFrmMtrx <- function(ndmtrx, all_ids, tids, spns, parallel) {
  res <- .getNdsPrdstsFrmMtrx(ndmtrx, all_ids, tids, spns,
    parallel = parallel, progress = "none"
  )
  max(res)
}

.newNd <- function(tree, id) {
  nd <- tree@ndlst[[id]]
  if (!tree@wspn) {
    spn <- pd <- prdst <- numeric()
  } else {
    spn <- nd[["spn"]]
    pd <- .getNdPDFrmLst(tree@ndlst,
      prinds = tree@prinds,
      id = id
    )
    prdst <- .getNdPrdstsFrmLst(tree@ndlst,
      prinds = tree@prinds,
      id = id
    )
  }
  if (is.null(nd[["txnym"]])) {
    txnym <- vector()
  } else {
    txnym <- nd[["txnym"]]
  }
  kids <- .getNdKidsFrmLst(tree@ndlst, prinds = tree@prinds, id = id, tinds = tree@tinds)
  new("Node",
    id = nd[["id"]], spn = spn, prid = as.character(nd[["prid"]][1]),
    ptid = as.character(nd[["ptid"]]), kids = as.character(kids),
    nkids = length(kids), pd = pd, txnym = txnym, prdst = prdst,
    root = tree@root == nd[["id"]], tip = length(nd[["ptid"]]) == 0
  )
}

#' @name Node-class
#' @aliases Node-method
#' @param x \code{Node} object
#' @param object \code{Node} object
#' @param i slot name
#' @param j missing
#' @param ... missing
#' @param drop missing
#' @title Node-class
#' @description The \code{Node} is an S4 class used for displaying node information.
#' It is only generated when a user implements the \code{[[]]} on a tree. Information
#' is only accurate if tree has been updated with \code{updateTree()}.
#' @slot id unique ID for node in tree['ndlst']
#' @slot spn length of preceding branch
#' @slot prid parent node ID
#' @slot ptid child node ID
#' @slot kids descending tip IDs
#' @slot nkids number of descending tip IDs
#' @slot txnym list of associated taxonyms
#' @slot pd total branch length represented by node
#' @slot prdst total branch length of connected prids
#' @slot root T/F root node?
#' @slot tip T/F tip node?
#' @exportClass Node
#' @seealso
#' \code{\link{TreeMan-class}}, \code{\link{TreeMen-class}}
setClass("Node",
  representation = representation(
    id = "character", # unique ID for node in tree@nodelist
    spn = "numeric", # length of preceding branch
    prid = "character", # parent node ID
    ptid = "vector", # child node IDs
    kids = "vector", # descending tip IDs
    nkids = "numeric", # number of descending tips
    txnym = "vector", # list of associated taxonyms
    pd = "numeric", # total branch length represented by node
    prdst = "numeric", # total branch length of connected pres
    root = "logical", # T/F root node?
    tip = "logical"
  ) # T/F tip node?
)
#' @rdname Node-class
#' @exportMethod as.character
setMethod(
  "as.character", c("x" = "Node"),
  function(x) {
    paste0("Node Obj. [ID=", x@id, "]")
  }
)
#' @rdname Node-class
#' @exportMethod show
setMethod(
  "show", "Node",
  function(object) {
    cat(summary(object))
  }
)
#' @rdname Node-class
#' @exportMethod print
setMethod(
  "print", "Node",
  function(x) {
    print(summary(x))
  }
)
#' @rdname Node-class
#' @exportMethod summary
setMethod(
  "summary", c("object" = "Node"),
  function(object) {
    if (object@root) {
      msg <- paste0("Node (root node):\n")
    } else if (object@tip) {
      msg <- paste0("Node (tip node):\n")
    } else {
      msg <- paste0("Node (internal node):\n")
    }
    msg <- paste0(msg, '  + ID: \"', object@id, '\"\n')
    if (length(object@txnym) > 0) {
      msg <- paste0(msg, '  + txnym: \"', paste0(object@txnym, collapse = '\", \"'), '\"\n')
    }
    if (!object@root) {
      msg <- paste0(msg, '  + prid: \"', object@prid, '\"\n')
    }
    if (!object@tip) {
      msg <- paste0(msg, '  + ptid: \"', paste0(object@ptid, collapse = '\", \"'), '\"\n')
      msg <- paste0(msg, "  + nkids: ", length(object@kids), "\n")
    }
    if (length(object@spn) > 0) {
      if (!object@root) {
        msg <- paste0(msg, "  + spn: ", signif(object@spn, 2), "\n")
      }
      msg <- paste0(msg, "  + predist: ", signif(object@prdst, 2), "\n")
      msg <- paste0(msg, "  + pd: ", signif(object@pd, 2), "\n")
    }
    cat(msg)
  }
)

#' @rdname Node-class
#' @exportMethod [
setMethod(
  "[", c("Node", "character", "missing", "missing"),
  function(x, i, j, ..., drop = TRUE) {
    if (!i %in% slotNames(x)) {
      stop(paste0(i, "  not in node"))
    }
    slot(x, i)
  }
)

#' @name writeTree
#' @title Write a Newick tree
#' @description Creates a Newick tree from a \code{TreeMan} object.
#' @details The \code{ndLabels} argument can be used to add a user defined node label in
#' the Newick tree. It should take only 1 argument, \code{nd}, the node represented as a list.
#' It should only return a single character value that can be added to a newick string.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @param append T/F append tree to already existing file
#' @param ndLabels node label function
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Newick_format},
#' \code{\link{readTree}}, \code{\link{randTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}},
#' \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # write out the tree with node labels as IDs
#' ndLabels <- function(n) {
#'   n[["id"]]
#' }
#' writeTree(tree, file = "example.tre", ndLabels = ndLabels)
#' file.remove("example.tre")
writeTree <- function(tree, file, append = FALSE, ndLabels = function(nd) {
                        return(NULL)
                      }, parallel = FALSE, progress = "none") {
  if (is(tree) == "TreeMen") {
    plyr::m_ply(tree@treelst,
      .fun = .writeTree, file = file,
      append = TRUE, ndLabels = ndLabels,
      .progress = progress, .parallel = parallel
    )
  } else if (is(tree) == "TreeMan") {
    .writeTree(tree, file, append, ndLabels)
  } else {
    stop("`tree` must be TreeMan or TreeMen")
  }
  NULL
}

.writeTree <- function(tree, file, append, ndLabels) {
  tipBytip <- function(i) {
    kids <- getNdKids(tree, prid)
    ids <- c(kids, prid, ndlst[[prid]][["prid"]])
    id <<- ids[!ids %in% deja_vues][1]
    deja_vues[i] <<- id
    spn <- ndlst[[id]][["spn"]]
    if (id %in% tids) {
      prids <- getNdPrids(tree, id)
      dpth <- which(prids == prid) - 1
      prid <<- ndlst[[id]][["prid"]]
      tpstr <- paste0(id, ":", spn)
      if (dpth > 0) {
        brckts <- paste0(rep("(", dpth), collapse = "")
        trstr <<- paste0(trstr, ",", brckts, tpstr)
      } else {
        trstr <<- paste0(trstr, ",", tpstr)
      }
    } else {
      prid <<- ndlst[[id]][["prid"]]
      ndlbl <- ndLabels(ndlst[[id]])
      trstr <<- paste0(trstr, ")", ndlbl, ":", spn)
    }
    NULL
  }
  # start with first tip
  # loop through tree structure adding tip by tip to string
  # unpack
  ndlst <- tree@ndlst
  tids <- tree@tips
  nids <- tree@nds
  rid <- tree@root
  # add first tip
  id <- tids[1]
  trstr <- ""
  deja_vues <- rep(NA, length(ndlst))
  deja_vues[1] <- id
  spn <- ndlst[[id]][["spn"]]
  dpth <- length(getNdPrids(tree, id))
  prid <- ndlst[[id]][["prid"]]
  tpstr <- paste0(id, ":", spn)
  trstr <- paste0(rep("(", dpth), collapse = "")
  trstr <- paste0(trstr, tpstr)
  # loop through nodes
  plyr::m_ply(2:(length(ndlst) - 1), .fun = tipBytip)
  ndlbl <- ndLabels(ndlst[[rid]])
  spn <- ndlst[[rid]][["spn"]]
  trstr <- paste0(trstr, ")", ndlbl, ":", spn, ";")
  write.table(
    x = trstr, file = file, quote = FALSE, row.names = FALSE,
    col.names = FALSE, append = append
  )
}

#' @name readTree
#' @title Read a Newick tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a Newick treefile
#' @details Read a single or multiple trees from a file, or a text string. Parallelizable
#' when reading multiple trees.
#' The function will add any internal node labels in the Newick tree as a user-defined data slots.
#' The name of this slot is defined with the \code{spcl_slt_nm}.
#' These data can be accessed/manipulated with the \code{`getNdsSlt()`} function.
#' Trees are always read as rooted. (Unrooted trees have polytomous root nodes.)
#' @param file file path
#' @param text Newick character string
#' @param spcl_slt_nm name of special slot for internal node labels, default 'Unknown'.
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Newick_format},
#' \code{\link{addNdmtrx}}, \code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{readTrmn}}, \code{\link{writeTrmn}},
#' \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#'
#' # tree string with internal node labels as bootstrap results
#' tree <- readTree(
#'   text = "((A:1.0,B:1.0)0.9:1.0,(C:1.0,D:1.0)0.8:1.0)0.7:1.0;",
#'   spcl_slt_nm = "bootstrap"
#' )
#' # retrieve bootstrap values by node
#' tree["bootstrap"]
readTree <- function(file = NULL, text = NULL, spcl_slt_nm = "Unknown", wndmtrx = FALSE,
                     parallel = FALSE, progress = "none") {
  if (!is.null(file)) {
    trstr <- scan(file, what = "raw", quiet = TRUE)
  } else {
    trstr <- text
  }
  if (length(trstr) > 1) {
    trstr <- as.list(trstr)
    trees <- plyr::mlply(trstr,
      .fun = .readTree, spcl_slt_nm = spcl_slt_nm,
      wndmtrx = wndmtrx, .progress = progress, .parallel = parallel
    )
    names(trees) <- NULL
    trees <- trees[1:length(trees)]
    tree <- as(trees, "TreeMen")
  } else {
    tree <- .readTree(trstr, spcl_slt_nm, wndmtrx)
  }
  tree
}

#' @useDynLib phylotaR
#' @useDynLib phylotaR cFindPrids
.readTree <- function(trstr, spcl_slt_nm, wndmtrx) {
  # Internals
  .idspn <- function(i) {
    mtdt <- substr(trstr, start = nds[i - 1] + 1, stop = nds[i])
    mtdt <- gsub("(\\(|\\)|,|;)", "", mtdt)
    mtdt <- strsplit(mtdt, ":")[[1]]
    id <- NA
    if (length(mtdt) == 0) {
      spn <- NA
    } else if (length(mtdt) == 1) {
      id <- mtdt
      spn <- NA
    } else if (length(mtdt) > 1 && mtdt[1] != "") {
      id <- mtdt[1]
      spn <- as.numeric(mtdt[2])
    } else {
      spn <- as.numeric(mtdt[2])
    }
    c(id, spn)
  }
  .add <- function(i) {
    nd <- vector("list", length = 4)
    names(nd) <- c("id", "ptid", "prid", "spn")
    nd[["id"]] <- ids[i]
    nd[["spn"]] <- spns[i]
    nd[["prid"]] <- ids[prinds[i]]
    nd[["ptid"]] <- ptids[ptnds_pool == i]
    nd
  }
  # get nodes from string
  nds <- c(1, as.integer(gregexpr("(,|\\))", trstr)[[1]]) - 1)
  nds <- c(nds, nchar(trstr))
  # get id and spn
  mtdt <- mapply(2:length(nds), FUN = .idspn)
  ids <- mtdt[1, ]
  spns <- as.numeric(mtdt[2, ])
  rm(mtdt)
  nds <- nds[-1]
  # gen prids
  opns <- gregexpr("\\(", trstr)[[1]]
  clss <- gregexpr("\\)", trstr)[[1]]
  prinds <- .Call("cFindPrids",
    PACKAGE = "phylotaR",
    as.integer(nds),
    as.integer(clss),
    as.integer(opns)
  )
  if (sum(prinds == -1) > 1) {
    stop("Invalid tree string")
  }
  root <- which(prinds == -1)
  prinds <- match(prinds, nds)
  tinds <- which(!1:length(ids) %in% prinds)
  prinds[is.na(prinds)] <- root
  spns[is.na(spns)] <- 0
  # move internal node labels to other
  other <- rep(NA, length(ids))
  intnds <- 1:length(ids) %in% prinds
  other[intnds] <- ids[intnds]
  # overwrite internal node ids only if at least one is malformed
  if (any(grepl("[^a-zA-Z_0-9]", ids[intnds]))) {
    ids[intnds] <- paste0("n", which(intnds))
  } else {
    # if using internal node ids as ids, no need for special slot
    other[intnds] <- NA
  }
  # rm NAs from IDs
  pull <- is.na(ids)
  ids[pull] <- paste0("n", which(pull))
  # ensure no dups in ids
  dups <- duplicated(ids)
  if (any(dups)) {
    dups <- unique(ids[dups])
    for (dup in dups) {
      pull <- ids == dup
      other[pull] <- ids[pull]
      ids[pull] <- paste0("n", which(pull))
    }
  }
  ptids <- ids[-root]
  ptnds_pool <- prinds[-root]
  ndlst <- lapply(1:length(ids), .add)
  names(ndlst) <- ids
  tree <- new("TreeMan",
    ndlst = ndlst, root = ids[root],
    ndmtrx = NULL, wtxnyms = FALSE,
    prinds = prinds, tinds = tinds
  )
  pull <- !is.na(other)
  if (any(pull)) {
    tree <- setNdsOther(tree,
      ids = ids[pull], vals = other[pull],
      slt_nm = spcl_slt_nm
    )
  }
  tree <- updateSlts(tree)
  if (wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name writeTrmn
#' @title Write a .trmn tree
#' @description Write to disk a \code{TreeMan} or \code{TreeMan} object using the .trmn treefile
#' @details Write a tree(s) to file using the .trmn format.
#' It is faster to read and write tree files using treeman with the .trmn file format.
#' In addition it is possible to encode more information than possible with the
#' Newick, e.g. any taxonomic information and additional slot names added to
#' the tree are recorded in the file.
#' @param tree TreeMan object or TreeMen object
#' @param file file path
#' @seealso
#' \code{\link{readTrmn}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' writeTrmn(tree, file = "test.trmn")
#' tree <- readTrmn("test.trmn")
#' file.remove("test.trmn")
writeTrmn <- function(tree, file) {
  .unpack <- function(ntree) {
    .makeDataFrame(ntree, tree@treelst[[ntree]])
  }
  .makeDataFrame <- function(ntree, tree) {
    res <- data.frame(tree = ntree, prind = tree@prinds)
    res[["id"]] <- names(tree@ndlst)
    if (tree@wspn) {
      res[["spn"]] <- vapply(
        tree@ndlst, function(x) x[["spn"]],
        numeric(1)
      )
    }
    if (tree@wtxnyms) {
      res[["txnym"]] <- vapply(
        tree@ndlst,
        function(x) paste0(x[["txnym"]], collapse = "|"),
        character(1)
      )
    }
    # add any additional slots
    if (length(tree@othr_slt_nms) > 0) {
      for (slt_nm in tree@othr_slt_nms) {
        res[[slt_nm]] <- sapply(
          tree@ndlst,
          function(x) x[[slt_nm]]
        )
      }
    }
    res
  }
  if ("TreeMan" %in% is(tree)) {
    res <- .makeDataFrame(1, tree)
  } else if ("TreeMen" %in% is(tree)) {
    res <- plyr::mdply(
      .data = data.frame(ntree = 1:tree@ntrees),
      .fun = .unpack
    )
    res <- res[, -1]
  } else {
    stop("`tree` must be TreeMan or TreeMen object.")
  }
  write.csv(res, file = file, quote = FALSE, row.names = FALSE)
}

#' @name readTrmn
#' @title Read a .trmn tree
#' @description Return a \code{TreeMan} or \code{TreeMen} object from a .trmn treefile
#' @details Read a tree(s) from a file using the .trmn format.
#' It is faster to read and write tree files using treeman with the .trmn file format.
#' In addition it is possible to encode more information than possible with the
#' Newick, e.g. any taxonomic information and additional slot names added to
#' the tree are recorded in the file.
#' @param file file path
#' @param wndmtrx T/F add node matrix? Default FALSE.
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{writeTrmn}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{randTree}}, \code{\link{saveTreeMan}}, \code{\link{loadTreeMan}}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' writeTrmn(tree, file = "test.trmn")
#' tree <- readTrmn("test.trmn")
#' file.remove("test.trmn")
readTrmn <- function(file, wndmtrx = FALSE, parallel = FALSE,
                     progress = "none") {
  .pack <- function(i) {
    .readTrmn(
      inpt[inpt[["tree"]] == i, ],
      wndmtrx
    )
  }
  inpt <- read.csv(file, stringsAsFactors = FALSE)
  trids <- unique(inpt[["tree"]])
  trees <- plyr::mlply(
    .data = trids, .fun = .pack,
    .parallel = parallel, .progress = progress
  )
  if (length(trees) == 1) {
    res <- trees[[1]]
  } else {
    trees <- trees[1:length(trees)]
    names(trees) <- NULL
    res <- as(trees, "TreeMen")
  }
  res
}

.readTrmn <- function(inpt, wndmtrx) {
  .add <- function(i) {
    nd <- vector("list", length = 4)
    names(nd) <- c("id", "ptid", "prid", "spn")
    nd[["id"]] <- ids[i]
    nd[["spn"]] <- spns[i]
    nd[["prid"]] <- ids[prinds[i]]
    nd[["ptid"]] <- ptids[ptnds_pool == i]
    nd
  }
  prinds <- inpt[["prind"]]
  # all internal nodes should occur more than once (twice for bifurcating trees)
  prind_test <- sum(prinds == 1:length(prinds)) == 1
  prind_test <- all(table(prinds) > 1) & prind_test
  if (!prind_test) {
    stop("Tree is corrupted, check node structure is hierarchical.")
  }
  ids <- inpt[["id"]]
  if ("spn" %in% names(inpt) && !is.na(inpt[["spn"]][1])) {
    spns <- inpt[["spn"]]
  } else {
    spns <- rep(0, length(ids))
  }
  tinds <- which(!1:length(ids) %in% prinds)
  root <- which(1:length(prinds) == prinds)
  ptids <- ids[-root]
  ptnds_pool <- prinds[-root]
  ndlst <- lapply(1:length(ids), .add)
  names(ndlst) <- ids
  tree <- new("TreeMan",
    ndlst = ndlst, root = ids[root],
    ndmtrx = NULL, wtxnyms = FALSE,
    prinds = prinds, tinds = tinds
  )
  if ("txnym" %in% names(inpt) && !is.na(inpt[["txnym"]][1])) {
    txnyms <- strsplit(inpt[["txnym"]], "\\|")
    names(txnyms) <- ids
    tree <- setTxnyms(tree, txnyms)
  }
  othr_slt_nms <- names(inpt)[!names(inpt) %in%
    c("id", "prind", "spn", "txnym", "tree")]
  if (length(othr_slt_nms) > 0) {
    for (slt_nm in othr_slt_nms) {
      tree <- setNdsOther(tree,
        ids = inpt[["id"]],
        vals = inpt[[slt_nm]], slt_nm = slt_nm
      )
    }
  }
  tree <- updateSlts(tree)
  if (wndmtrx) {
    tree <- addNdmtrx(tree)
  }
  tree
}

#' @name saveTreeMan
#' @title Save a TreeMan object in serialization format
#' @description \code{TreeMan} equivalent to \code{save()} but able to handle
#' node matrices.
#' @details It is not possible to use \code{save()} on \code{TreeMan} objects
#' with node matrices. Node matrices are bigmemory matrices and are therefore outside
#' the R environment, see bigmemory documentation for more information. Saving and loading
#' a bigmemory matrix may cause memory issues in R and cause R to crash.
#'
#' This function can safely store a \code{TreeMan} object with and without
#' a node matrix. This function stores the tree using the serialization format and the node
#' matrix as a hidden .csv. Both parts of the tree can be reloaded to an R environment
#' with \code{loadTreeMan()}. The hidden node matrix filename is based on the file argument:
#' \code{file + _ndmtrx}
#'
#' Reading and writing trees with \code{saveTreeMan()} and
#' \code{loadTreeMan} is faster than any of the other read and write functions.
#' @param tree \code{TreeMan} object
#' @param file file path
#' @seealso
#' \code{\link{loadTreeMan}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}}
#' @export
#' @examples
#'
#' tree <- randTree(100, wndmtrx = TRUE)
#' saveTreeMan(tree, file = "test.RData")
#' rm(tree)
#' tree <- loadTreeMan(file = "test.RData")
#' file.remove("test.RData", "testRData_ndmtrx")
saveTreeMan <- function(tree, file) {
  ndmtrx_file <- paste0(gsub("\\.", "", file), "_ndmtrx")
  if (!is.null(tree@ndmtrx)) {
    bigmemory::write.big.matrix(x = tree@ndmtrx, filename = ndmtrx_file)
    tree <- rmNdmtrx(tree)
  }
  save(list = c("tree", "ndmtrx_file"), file = file)
}

#' @name loadTreeMan
#' @title Load a TreeMan object in serialization format
#' @description \code{TreeMan} equivalent to \code{load()} but able to handle
#' node matrices.
#' @details It is not possible to use \code{save()} on \code{TreeMan} objects
#' with node matrices. Node matrices are bigmemory matrices and are therefore outside
#' the R environment, see bigmemory documentation for more information. Saving and loading
#' a bigmemory matrix may cause memory issues in R and cause R to crash.
#'
#' This function can safely read a \code{TreeMan} object with and without
#' a node matrix. \code{saveTreeMan()} function stores the tree using the serialization format
#' and the node matrix as a hidden .csv. Both parts of the tree can be reloaded to an R environment
#' with \code{loadTreeMan()}. The hidden node matrix filename is based on the file argument:
#' \code{file + _ndmtrx}
#'
#' Reading and writing trees with \code{saveTreeMan()} and
#' \code{loadTreeMan} is faster than any of the other read and write functions.
#' @param file file path
#' @seealso
#' \code{\link{saveTreeMan}},
#' \code{\link{readTree}},\code{\link{writeTree}},
#' \code{\link{readTrmn}}, \code{\link{writeTrmn}}
#' @export
#' @examples
#'
#' tree <- randTree(100, wndmtrx = TRUE)
#' saveTreeMan(tree, file = "test.RData")
#' rm(tree)
#' tree <- loadTreeMan(file = "test.RData")
#' file.remove("test.RData", "testRData_ndmtrx")
loadTreeMan <- function(file) {
  ndmtrx_file <- NULL
  load(file)
  if (file.exists(ndmtrx_file)) {
    tree@ndmtrx <- bigmemory::read.big.matrix(
      filename = ndmtrx_file,
      type = "integer", shared = FALSE
    )
  }
  tree
}

#' @name searchTxnyms
#' @title Get node labels based on online taxonomic database
#' @description Return names of each node in tree based on searching tip labels
#'  through Global Names Resolver \url{http://resolver.globalnames.org/} in NCBI.
#' @details For each node, all the descendants are searched, the taxonomic lineages returned and
#' then searched to find the lowest shared name.
#' All the tip labels are searched against a specified taxonomic database through the GNR and NCBI.
#' (So far only tested with NCBI database.)
#' Use the infer argument to ensure a taxonym is returned for all nodes. If infer is true,
#' all nodes without an identifed taxonym adopt the taxonym of their parent.
#' Will raise a warning if connection fails and will return NULL.
#' @param tree TreeMan object
#' @param cache T/F, create a local cache of downloaded names?
#' @param parent specify parent of all names to prevent false names
#' @param clean T/F, ensure returned names contain no special characters?
#' @param infer T/F, infer taxonyms for unfound nodes?
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{setTxnyms}}, \code{\link{getNdsFrmTxnyms}}
#' @export
#' @examples
#' tree <- randTree(8)
#' new_tids <- c(
#'   "Gallus_gallus", "Aileuropoda_melanoleucha", "Ailurus_fulgens",
#'   "Rattus_rattus", "Mus_musculus", "Gorilla_gorilla", "Pan_trogoldytes", "Homo_sapiens"
#' )
#' tree <- setNdsID(tree, tree["tips"], new_tids)
#' nd_labels <- searchTxnyms(tree)
#' print(nd_labels)
# TODO: add compatibility with other GNR datasources
# TODO: catalogue of life, unlike NCBI, does not keep lineages and rank lengths constant between names
searchTxnyms <- function(tree, cache = FALSE, parent = NULL, clean = TRUE,
                         infer = TRUE) {
  # Use GNR to label all nodes in a phylogeny
  # first replace all _ with spaces
  tip_labels <- gsub("_", " ", tree@tips)
  nids <- tree@nds
  nd_labels <- rep(NA, tree@nall)
  names(nd_labels) <- tree@all
  taxa_res <- taxaResolve(tip_labels, datasource = 4, cache = cache, parent = parent)
  if (is.null(taxa_res)) {
    return(NULL)
  }
  # for tips use the first word of the name
  nd_labels[tree@tips] <- vapply(
    strsplit(tip_labels, "\\s+"), function(x) x[1],
    character(1)
  )
  nds_kids <- getNdsKids(tree, ids = nids)
  for (nid in nids) {
    kids <- gsub("_", " ", nds_kids[[nid]])
    genus_names <- vapply(
      strsplit(kids, "\\s+"), function(x) x[1],
      character(1)
    )
    if (all(genus_names == genus_names[1])) {
      nd_labels[nid] <- genus_names[1]
    } else {
      lineages <- as.character(taxa_res[taxa_res$search.name %in% kids, "lineage"])
      lineages <- strsplit(lineages, "\\|")
      lineages <- lineages[!is.na(lineages)]
      if (length(lineages) > 1) {
        nd_labels[nid] <- .findClade(lineages)
      }
    }
  }
  if (clean) {
    nd_labels <- gsub("\\s", "_", nd_labels)
    nd_labels <- gsub("[^a-zA-Z_0-9]", "", nd_labels)
  }
  if (infer & any(is.na(nd_labels))) {
    if (sum(is.na(nd_labels)) / length(nd_labels) > 0.5) {
      message(
        "Fewer than 50% of nodes identified,",
        " not attempting inference of remaning nodes."
      )
    } else {
      for (i in which(is.na(nd_labels))) {
        prids <- getNdPrids(tree, names(nd_labels)[i])
        pssbls <- nd_labels[prids]
        pssbls <- pssbls[!is.na(pssbls)]
        if (length(pssbls) > 0) {
          nd_labels[[i]] <- pssbls[[1]]
        }
      }
    }
  }
  nd_labels
}

#' @name taxaResolve
#' @title Resolve taxonmic names online
#' @description Resolve taxonomic names via the Global Names Resolver.
#' @details Returns dataframe containing GNR metadata for each name wames
#' that cannot be resolved are returned as NA. Various datasources are
#' available, see \url{http://resolver.globalnames.org/data_sources} for a
#' list and IDs. Default is 4 for NCBI.
#' Will raise a warning if connection fails and will return NULL.
#' @param nms vector of names
#' @param batch size of the batches to be queried
#' @param datasource ID number of the datasource
#' @param genus boolean, if true will search against GNR with just the genus
#'  name for names that failed to resolve using the full species name
#' @param cache T/F, create a local cache of downloaded names?
#' @param parent specify parent of all names to prevent false names
#' @seealso
#' \code{\link{searchTxnyms}}, \code{\link{setTxnyms}}, \code{\link{getNdsFrmTxnyms}}
#' @export
#' @examples
#' my_lovely_names <- c(
#'   "Gallus gallus", "Pongo pingu", "Homo sapiens",
#'   "Arabidopsis thaliana", "Macaca thibetana", "Bacillus subtilis"
#' )
#' res <- taxaResolve(nms = my_lovely_names)
#' length(colnames(res)) # 10 different metadata for returned names including original search name
#' # let's look at the lineages
#' lineages <- strsplit(as.vector(res$lineage), "\\|")
#' print(lineages[[6]]) # the bacteria has far fewer taxonomic levels
# NOTE. Originally built for MTT
taxaResolve <- function(nms, batch = 100, datasource = 4, genus = TRUE,
                        cache = FALSE, parent = NULL) {
  .replace <- function(i, slot.name) {
    # controlled extraction
    element <- data[[i]]$result[[1]][[slot.name]]
    if (!is.null(element)) {
      res <- element
    } else {
      res <- NA
    }
    res
  }
  batchResolve <- function(batch.nms) {
    # create query from nms
    url <- "http://resolver.globalnames.org/name_resolvers.json?"
    data_source_ids <- paste0("&data_source_ids=", datasource)
    nms2 <- paste0("names=", paste0(stringr::str_replace_all(
      batch.nms, " ", "+"
    ), collapse = "|"))
    query <- paste(plyr::compact(list(
      url, nms2,
      data_source_ids
    )), collapse = "")
    # search via API
    data <- .safeFromJSON(query)$data
    return(data)
  }
  # avoid names -- names exist on database but mean nothing
  avoid <- c("unidentified")
  # make sure names don't have '_'
  trms <- gsub("_", " ", nms)
  # remove all non alphanumerics
  trms <- gsub("\\W+", " ", trms)
  # remove trailing whitespace
  trms <- gsub("^\\s+|\\s+$", "", trms)
  # any missing trms, replace with stubs
  trms[trms == ""] <- "invalid"
  deja_vues <- rep(FALSE, length(trms))
  data <- vector("list", length = length(trms))
  names(data) <- nms
  if (cache) {
    if (!file.exists("gnr_cache")) {
      dir.create("gnr_cache")
    }
    for (i in 1:length(nms)) {
      fp <- file.path("gnr_cache", paste0(nms[i], ".RData"))
      if (file.exists(fp)) {
        load(fp)
        data[[nms[i]]] <- nd
        deja_vues[i] <- TRUE
      }
    }
  }
  if (sum(!deja_vues) > 0) {
    # Split nms into batch sized chunks
    #  http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
    x <- seq_along(trms[!deja_vues])
    btrms <- split(trms[!deja_vues], ceiling(x / batch))
    bnms <- split(nms[!deja_vues], ceiling(x / batch))
    for (i in 1:length(btrms)) {
      temp.data <- batchResolve(btrms[[i]])
      if (is.null(temp.data)) {
        return(NULL)
      }
      data[bnms[[i]]] <- temp.data
    }
  }
  # transform results into output
  search.name <- name.string <- canonical.form <-
    lineage <- lineage.ids <- rank <- taxid <-
    match.type <- prescore <- score <- rep(NA, length(nms))
  for (i in 1:length(data)) {
    parent_test <- TRUE
    nd <- data[[i]]
    if (cache & !deja_vues[i]) {
      fp <- file.path("gnr_cache", paste0(names(data)[i], ".RData"))
      save(nd, file = fp)
    }
    if (!"results" %in% names(nd)) {
      search.name[i] <- nms[i]
    } else if (nd[[1]] %in% avoid) {
      search.name[i] <- nms[i]
    } else {
      search.name[i] <- nms[i]
      lng <- .replace(i, "classification_path")
      if (!is.null(parent)) {
        parent_test <- grepl(parent, lng)
      }
      if (parent_test) {
        name.string[i] <- .replace(i, "name_string")
        canonical.form[i] <- .replace(i, "canonical_form")
        lineage[i] <- lng
        lineage.ids[i] <- .replace(i, "classification_path_ids")
        rank[i] <- .replace(i, "classification_path_ranks")
        taxid[i] <- .replace(i, "taxon_id")
        match.type[i] <- .replace(i, "match_type")
        prescore[i] <- .replace(i, "prescore")
        score[i] <- nd$results[[1]]$score
      }
    }
  }
  res <- data.frame(
    search.name = search.name,
    name.string = name.string,
    canonical.form = canonical.form,
    lineage = lineage, lineage.ids = lineage.ids,
    rank = rank, taxid = taxid,
    match.type = match.type, prescore = prescore,
    score = score, stringsAsFactors = FALSE
  )
  failed <- which(is.na(res$name.string))
  if (genus & length(failed) > 0) {
    # if genus, search just genus names
    genus.nms <- sub("\\s+.*", "", res$search.name[failed])
    genus.res <- taxaResolve(genus.nms, batch, datasource,
      genus = FALSE,
      parent = parent, cache = cache
    )
    # replace in original results, all slots except search.name
    res[failed, -1] <- genus.res[, -1]
  }
  return(res)
}

.safeFromJSON <- function(url, max_trys = 5, power = 2) {
  # Safe wrapper for fromJSON
  trys <- 0
  waittime <- 2
  while (trys < max_trys) {
    json_obj <- try(RJSONIO::fromJSON(url), silent = TRUE)
    if (inherits(json_obj, "try-error") {
      cat("---- Connection failed: trying again in [", waittime,
        "s]----\n",
        sep = ""
      )
      trys <- trys + 1
      Sys.sleep(waittime)
      waittime <- waittime * power
    } else {
      return(json_obj)
    }
  }
  warning("Failed to connect, server may be down.")
  list("data" = NULL)
}

.findClade <- function(lineages) {
  # for a list of lineages, find the clade shared by all
  subj <- lineages[[1]]
  for (i in 2:length(lineages)) {
    query <- lineages[[i]]
    subj <- subj[subj %in% query]
  }
  subj[length(subj)]
}

#' @name setTxnyms
#' @title Set the txnym slots in a tree
#' @description Return a tree with txnyms added to specified nodes
#' @details Returns a tree. Specify the taxonomic groups for nodes in a tree
#' by providing a vector or list named by node IDs. Takes output from \code{searchTxnyms}.
#' Only letters, numbers and underscores allowed. To remove special characters use regular
#' expressions, e.g. \code{gsub(['a-zA-Z0-9_'], '', txnym)}
#' @param tree \code{TreeMan} object
#' @param txnyms named vector or list
#' @seealso
#' \code{\link{taxaResolve}}, \code{\link{searchTxnyms}},
#' \code{\link{getNdsLng}}, \code{\link{getNdLng}},
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' data(mammals)
#' # let's change the txnym for humans
#' # what's its summary before we change anything?
#' summary(mammals[["Homo_sapiens"]])
#' # now let's add Hominini
#' new_txnym <- list("Homo_sapiens" = c("Hominini", "Homo"))
#' mammals <- setTxnyms(mammals, new_txnym)
#' summary(mammals[["Homo_sapiens"]])
setTxnyms <- function(tree, txnyms) {
  .add <- function(nid) {
    for (txnym in txnyms[[nid]]) {
      if (grepl("[^a-zA-Z_0-9]", txnym)) {
        stop(paste0(
          "Unsuitable characters in [",
          txnym, "]"
        ))
      }
    }
    tree@ndlst[[nid]][["txnym"]] <<- txnyms[[nid]]
  }
  pull <- names(txnyms) %in% names(tree@ndlst)
  txnyms <- txnyms[pull]
  plyr::m_ply(names(txnyms), .fun = .add)
  tree@wtxnyms <- TRUE
  tree
}

#' @name setPD
#' @title Set the phylogenetic diversity
#' @description Return a tree with the phylogenetic diversity altered.
#' @details Use this function to convert the phylogenetic diversity of a tree. For example,
#' you might want to convert the tree so the sum of all branches is 1. This function will achieve
#' that by modiyfing every branch, while maintaining their relative lengths.
#' @param tree \code{TreeMan} object
#' @param val new phylogenetic diversity
#' @seealso
#' \code{\link{setAge}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- setPD(tree, val = 1)
#' summary(tree)
setPD <- function(tree, val) {
  spns <- getNdsSlt(tree, ids = tree@all, slt_nm = "spn")
  spns <- spns / (tree@pd / val)
  tree <- setNdsSpn(tree, ids = tree@all, vals = spns)
  tree@pd <- val
  tree
}

#' @name setAge
#' @title Set the age of a tree
#' @description Return a tree with the age altered.
#' @details Use this function to change the age of a tree. For example,
#' you might want to convert the tree so that its age equals 1. This function will achieve
#' that by modiyfing every branch, while maintaining their relative lengths.
#' @param tree \code{TreeMan} object
#' @param val new age
#' @seealso
#' \code{\link{setPD}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- setAge(tree, val = 1)
#' summary(tree)
setAge <- function(tree, val) {
  # cdmuir's correction
  tree_age <- getAge(tree)
  spns <- getNdsSlt(tree, ids = tree@all, slt_nm = "spn")
  # spns <- spns/(tree_age) # previous, incorrect
  spns <- spns * val / tree_age
  tree <- setNdsSpn(tree, ids = tree@all, vals = spns)
  tree
}

#' @name setNdSpn
#' @title Set the branch length of a specific node
#' @description Return a tree with the span of a node altered.
#' @details Takes a tree, a node ID and a new value for the node's preceding branch length (span).
#' @param tree \code{TreeMan} object
#' @param id id of node whose preceding edge is to be changed
#' @param val new span
#' @seealso
#' \code{\link{setNdsSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- setNdSpn(tree, id = "t1", val = 100)
#' tree <- updateSlts(tree)
#' summary(tree)
setNdSpn <- function(tree, id, val) {
  tree@ndlst[[id]][["spn"]] <- val
  tree@updtd <- FALSE
  tree
}

#' @name setNdsSpn
#' @title Set the branch lengths of specific nodes
#' @description Return a tree with the spans of nodes altered.
#' @details Runs \code{setNdSpn} over multiple nodes. Parallelizable.
#' @param tree \code{TreeMan} object
#' @param ids ids of nodes whose preceding edges are to be changed
#' @param vals new spans
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdSpn}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # make tree taxonomic
#' tree <- setNdsSpn(tree, ids = tree["all"], vals = 1)
#' summary(tree)
#' # remove spns by setting all to 0
#' tree <- setNdsSpn(tree, ids = tree["all"], vals = 0)
#' summary(tree)
setNdsSpn <- function(tree, ids, vals, parallel = FALSE, progress = "none") {
  .reset <- function(id, spn) {
    ndlst[[id]][["spn"]] <- spn
    ndlst[[id]]
  }
  ndlst <- tree@ndlst[ids]
  l_data <- data.frame(id = ids, spn = vals, stringsAsFactors = FALSE)
  ndlst <- plyr::mlply(l_data,
    .fun = .reset, .parallel = parallel,
    .progress = progress
  )
  ndlst <- ndlst[1:length(ndlst)]
  tree@ndlst[ids] <- ndlst
  tree <- updateSlts(tree)
  tree
}

#' @name setNdID
#' @title Set the ID of a node
#' @description Return a tree with the ID of a node altered.
#' @details IDs cannot be changed directly for the \code{TreeMan} class. To change an
#' ID use this function. Warning: all IDs must be unique, avoid spaces in IDs and only
#' use letters, numbers and underscores.
#' Use \code{\link{updateSlts}} after running.
#' @param tree \code{TreeMan} object
#' @param id id to be changed
#' @param val new id
#' @seealso
#' \code{\link{setNdsID}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- setNdID(tree, "t1", "heffalump")
#' tree <- updateSlts(tree)
setNdID <- function(tree, id, val) {
  tree@updtd <- FALSE
  setNdsID(tree, id, val)
}

#' @name setNdsID
#' @title Set the IDs of multiple nodes
#' @description Return a tree with the IDs of nodes altered.
#' @details Runs \code{setNdID()} over multiple nodes. Warning: all IDs must be unique,
#' avoid spaces in IDs, only use numbers, letters and underscores. Parellizable.
#' @param tree \code{TreeMan} object
#' @param ids ids to be changed
#' @param vals new ids
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdID}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' new_ids <- paste0("heffalump_", 1:tree["ntips"])
#' tree <- setNdsID(tree, tree["tips"], new_ids)
#' summary(tree)
setNdsID <- function(tree, ids, vals, parallel = FALSE, progress = "none") {
  # internals
  .testSpcls <- function(id) {
    if (grepl("[^a-zA-Z_0-9]", id)) {
      stop(paste0("Unsuitable characters in [", id, "]"))
    }
    NULL
  }
  .rplcS4 <- function(slt) {
    if (any(slot(tree, slt) %in% ids)) {
      mtchs <- match(slot(tree, slt), ids)
      return(vals[mtchs])
    } else {
      return(slot(tree, slt))
    }
  }
  .reset <- function(i) {
    .rplc <- function(slt) {
      res <- nd[[slt]]
      mtchs <- match(res, ids)
      res[which(!is.na(mtchs))] <-
        vals[mtchs[!is.na(mtchs)]]
      res
    }
    nd <- tree@ndlst[[i]]
    nd[["id"]] <- .rplc("id")
    nd[["ptid"]] <- .rplc("ptid")
    nd[["prid"]] <- .rplc("prid")
    nd
  }
  mapply(FUN = .testSpcls, vals)
  l_data <- data.frame(i = 1:length(tree@ndlst), stringsAsFactors = FALSE)
  ndlst <- plyr::mlply(l_data, .fun = .reset, .parallel = parallel, .progress = progress)
  ndlst <- ndlst[1:length(ndlst)]
  all <- names(tree@ndlst)
  all[match(ids, all)] <- vals
  names(ndlst) <- all
  tree@ndlst <- ndlst
  tree@tips <- .rplcS4("tips")
  tree@nds <- .rplcS4("nds")
  tree@root <- .rplcS4("root")
  tree <- updateSlts(tree)
  tree
}

#' @name setNdOther
#' @title Set a user defined slot
#' @description Return a tree with a user defined slot for node ID.
#' @details A user can specify new slots in a tree. Add a new slot with this function
#' by providing a node ID, a value for the new slot and a unique new slot name. Slot names
#' must not be default \code{TreeMan} names. The new value can be any data type.
#' @param tree \code{TreeMan} object
#' @param id id of the node
#' @param val data for slot
#' @param slt_nm slot name
#' @seealso
#' \code{\link{setNdsOther}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' tree <- setNdOther(tree, "t1", 1, "binary_val")
#' tree <- updateSlts(tree)
#' (getNdSlt(tree, id = "t1", slt_nm = "binary_val"))
setNdOther <- function(tree, id, val, slt_nm) {
  tree@ndlst[[id]][slt_nm] <- val
  tree@updtd <- FALSE
  if (!slt_nm %in% tree@othr_slt_nms) {
    tree@othr_slt_nms <- c(tree@othr_slt_nms, slt_nm)
  }
  tree
}

#' @name setNdsOther
#' @title Set a user defined slot for multiple nodes
#' @description Return a tree with a user defined slot for node IDs.
#' @details Runs \code{setNdOther()} over multiple nodes. Parellizable.
#' @param tree \code{TreeMan} object
#' @param ids id sof the nodes
#' @param vals data for slot
#' @param slt_nm slot name
#' @param parallel logical, make parallel?
#' @param progress name of the progress bar to use, see \code{\link{create_progress_bar}}
#' @seealso
#' \code{\link{setNdOther}}
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' # e.g. confidences for nodes
#' vals <- runif(min = 0, max = 1, n = tree["nall"])
#' tree <- setNdsOther(tree, tree["all"], vals, "confidence")
#' tree <- updateSlts(tree)
#' summary(tree)
#' (getNdsSlt(tree, ids = tree["all"], slt_nm = "confidence"))
setNdsOther <- function(tree, ids, vals, slt_nm, parallel = FALSE, progress = "none") {
  .set <- function(id, val) {
    tree@ndlst[[id]][[slt_nm]] <<- val
  }
  l_data <- data.frame(
    id = ids, val = vals,
    stringsAsFactors = FALSE
  )
  plyr::m_ply(.data = l_data, .fun = .set, .parallel = parallel, .progress = progress)
  tree@updtd <- FALSE
  if (!slt_nm %in% tree@othr_slt_nms) {
    tree@othr_slt_nms <- c(tree@othr_slt_nms, slt_nm)
  }
  tree
}

#' @name rmOtherSlt
#' @title Remove a user-defined slot
#' @description Returns a tree with a user-defined tree slot removed.
#' @details A user can specify a new slot using the \code{setNdSlt()} function
#' or upon reading a tree. This can be removed using this function by specifying
#' the name of the slot to be removed.
#' @param tree \code{TreeMan} object
#' @param slt_nm name of slot to be removed
#' @seealso
#' \code{\link{setNdOther}}, \code{\link{setNdsOther}},
#' \url{https://github.com/DomBennett/treeman/wiki/set-methods}
#' @export
#' @examples
#'
#' tree <- randTree(10)
#' vals <- runif(min = 0, max = 1, n = tree["nall"])
#' tree <- setNdsOther(tree, tree["all"], vals, "confidence")
#' tree <- updateSlts(tree)
#' summary(tree)
#' tree <- rmOtherSlt(tree, "confidence")
#' tree <- updateSlts(tree)
#' summary(tree)
rmOtherSlt <- function(tree, slt_nm) {
  .set <- function(id) {
    tree@ndlst[[id]][[slt_nm]] <<- NULL
  }
  l_data <- data.frame(id = tree@all, stringsAsFactors = FALSE)
  plyr::m_ply(.data = l_data, .fun = .set)
  tree@updtd <- FALSE
  tree@othr_slt_nms <- tree@othr_slt_nms[tree@othr_slt_nms != slt_nm]
  tree
}

# roxygen imports
#' @import methods
#' @importFrom graphics lines plot.default text
#' @importFrom utils combn write.table read.csv write.csv
#' @importFrom stats runif

#' @name TreeMan-class
#' @title TreeMan-class
#' @aliases TreeMan-method
#' @description S4 class for representing phylogenetic trees as a list of nodes.
#' @param x \code{TreeMan} object
#' @param i node ID or slot name
#' @param object \code{TreeMan} object
#' @param max.level \code{str()} maximum number of levels to show
#' @param ... additional tree objects
#' @param j missing
#' @param drop missing
#' @slot ndlst list of nodes
#' @slot nds vector of node ids that are internal nodes
#' @slot nnds numeric of number of internal nodes in tree
#' @slot tips vector of node ids that are tips
#' @slot ntips numeric of number of internal nodes in tree
#' @slot all vector of all node ids
#' @slot nall numeric of number of all nodes in tree
#' @slot pd numeric of total branch length of tree
#' @slot tinds indexes of all tip nodes in tree
#' @slot prinds indexes of all pre-nodes in tree
#' @slot wspn logical, do nodes have spans
#' @slot wtxnyms logical, do nodes have txnyms
#' @slot ply logical, is tree bifurcating
#' @slot root character of node id of root, if no root then empty character
#' @slot updtd logical, if tree slots have been updated since initiation or change
#' @slot othr_slt_nms vector, character list of additional data slots added to nodes
#' @slot ndmtrx matrix, T/Fs representing tree structure
#' @details
#' A \code{TreeMan} object holds a list of nodes. The idea of the \code{TreeMan}
#' class is to make adding and removing nodes as similar as possible to adding
#' and removing elements in a list. Note that internal nodes and tips are
#' both considered nodes. Trees can be polytomous but not unrooted.
#'
#'
#' Each node within the \code{TreeMan} \code{ndlst} contains the following data slots:
#' \itemize{
#'    \item \code{id}: character string for the node ID
#'    \item \code{txnym}: name of taxonomic clade (optional)
#'    \item \code{spn}: length of the preceding branch
#'    \item \code{prid}: ID of the immediately preceding node, NULL if root
#'    \item \code{ptid}: IDs of the immediately connecting nodes
#' }
#'
#' See below in 'Examples' for these methods in use.
#' @seealso
#' \code{\link{randTree}}, \code{\link{Node-class}},
#' \code{\link{phylo-to-TreeMan}}, \code{\link{TreeMan-to-phylo}}
#' @examples
#'
#' # Generate random tree
#' tree <- randTree(10)
#' # Print to get basic stats
#' summary(tree)
#' # Slots....
#' tree["tips"] # return all tips IDs
#' tree["nds"] # return all internal node IDs
#' tree["ntips"] # count all tips
#' tree["nnds"] # count all internal nodes
#' tree["root"] # identify root node
#' tree[["t1"]] # return t1 node object
#' tree["pd"] # return phylogenetic diversity
#' tree["ply"] # is polytomous?
#' # Additional special slots (calculated upon call)
#' tree["age"] # get tree's age
#' tree["ultr"] # determine if tree is ultrametric
#' tree["spns"] # get all the spans of the tree IDs
#' tree["prids"] # get all the IDs of preceding nodes
#' tree["ptids"] # get all the IDs of following nodes
#' tree["txnyms"] # get all the taxonyms of all nodes
#' # In addition [] can be used for any user-defined slot
#' # Because all nodes are lists with metadata we can readily
#' #  get specific information on nodes of interest
#' nd <- tree[["n2"]]
#' summary(nd)
#' # And then use the same syntax for the tree
#' nd["nkids"] # .... nkids, pd, etc.
#'
#' # Convert to phylo and plot
#' library(ape)
#' tree <- as(tree, "phylo")
#' plot(tree)
#' @exportClass TreeMan
setClass("TreeMan",
  representation = representation(
    ndlst = "list", # list of node lists
    nds = "vector", # vector of node ids that are internal nodes
    nnds = "numeric", # numeric of number of internal nodes in tree
    tips = "vector", # vector of node ids that are tips
    ntips = "numeric", # numeric of number of internal nodes in tree
    all = "vector", # vector of all Node ids
    nall = "numeric", # numeric of number of all nodes in tree
    pd = "numeric", # numeric of total branch length of tree
    wspn = "logical", # logical, do all nodes have spans
    wtxnyms = "logical", # logical, do nodes txnyms
    ply = "logical", # logical, is tree bifurcating
    updtd = "logical", # logical, if tree slots has been updated since a change
    ndmtrx = "ANY", # bigmemory matrix of logicals
    tinds = "vector", # indexes of tip nodes
    prinds = "vector", # indexes of pre-nodes
    root = "character", # character of node id of root, if no root then empty character
    othr_slt_nms = "vector"
  ), # if new slots added to node, list them here
  validity = fastCheckTreeMan
)

# Accessor methods
#' @rdname TreeMan-class
#' @exportMethod [[
setMethod(
  "[[", c("TreeMan", "character"),
  function(x, i) {
    if (!i %in% names(x@ndlst)) {
      srch_trm <- gsub(" ", "_", i) # usual mistake
      pssbls <- which(agrepl(srch_trm, names(x@ndlst),
        ignore.case = TRUE,
        max.distance = 0.25
      ))
      pssbls <- names(x@ndlst)[pssbls]
      if (length(pssbls) > 0 & length(pssbls) < 50) {
        msg <- paste0("Can't find [", i, "]. Did you mean ....\n")
        for (p in pssbls) {
          msg <- paste0(msg, '"', p, '"\n')
        }
        msg <- paste0(msg, "?\n")
      } else {
        msg <- paste0("Can't find [", i, "] in tree.")
      }
      stop(msg)
    }
    .newNd(x, i)
  }
)
#' @rdname TreeMan-class
setMethod(
  "[", c("TreeMan", "character", "missing", "missing"),
  function(x, i, j, ..., drop = TRUE) {
    slt_nms <- slotNames(x)
    slt_nms <- slt_nms[slt_nms != "ndlst"]
    slt_nms <- slt_nms[slt_nms != "ndmtrx"]
    slt_nms <- slt_nms[slt_nms != "tinds"]
    slt_nms <- slt_nms[slt_nms != "prinds"]
    # ultr is special, shouldn't be updated when updateSlts()
    # too slow to calculate. Instead only calc if called.
    if (i == "ultr") {
      return(isUltrmtrc(x))
    }
    if (i == "age") {
      return(getAge(x))
    }
    # getNdsSlt extractor
    xtr_slts <- c("spns", "prids", "ptids", "txnyms")
    if (i %in% xtr_slts) {
      slt_nm <- sub("s$", "", i) # rm s at end
      res <- getNdsSlt(x, slt_nm, x@all)
      names(res) <- x@all
      return(res)
    }
    if (i %in% x@othr_slt_nms) {
      res <- getNdsSlt(x, i, x@all)
      names(res) <- x@all
      return(res)
    }
    if (!i %in% slt_nms) {
      slt_nms <- paste0(c(
        slt_nms, "ultr", "age", xtr_slts,
        x@othr_slt_nms
      ), collapse = ", ")
      stop(paste0(
        "`", i, "` not a tree slot. Available slots: ",
        slt_nms
      ))
    }
    slot(x, i)
  }
)

# display methods
#' @rdname TreeMan-class
#' @exportMethod as.character
setMethod(
  "as.character", c("x" = "TreeMan"),
  function(x) {
    paste0("TreeMan Object of [", length(x@tips), "] tips")
  }
)
#' @rdname TreeMan-class
#' @exportMethod show
setMethod(
  "show", "TreeMan",
  function(object) {
    msg <- as.character(object)
    cat(msg)
  }
)
#' @rdname TreeMan-class
#' @exportMethod print
setMethod(
  "print", "TreeMan",
  function(x) {
    msg <- as.character(x)
    print(msg)
  }
)
#' @rdname TreeMan-class
#' @exportMethod str
setMethod(
  "str", c("object" = "TreeMan"),
  function(object, max.level = 2L, ...) {
    if (is.na(max.level)) {
      stop("max.level must be numeric")
    }
    str@default(object, max.level = max.level, ...)
  }
)
#' @rdname TreeMan-class
#' @exportMethod summary
setMethod(
  "summary", c("object" = "TreeMan"),
  function(object) {
    if (!fastCheckTreeMan(object)) {
      stop("Tree is corrupted. Run `checkNdlst()` to see how.")
    }
    if (!object@updtd) {
      stop("Tree is not updated since change or initiation. Use `updateSlts()`")
    }
    msg <- "Tree (TreeMan Object):\n"
    msg <- paste0(msg, "  + ", object@ntips, " tips\n")
    msg <- paste0(msg, "  + ", object@nnds, " internal nodes\n")
    if (!is.null(object@ndmtrx)) {
      msg <- paste0(msg, "  + With node matrix\n")
    }
    if (object@wtxnyms) {
      msg <- paste0(msg, "  + With taxonomic names\n")
    }
    if (object@ply) {
      msg <- paste0(msg, "  + Polytomous\n")
    } else {
      msg <- paste0(msg, "  + Binary\n")
    }
    if (length(object@root) == 0) {
      if (!object@wspn) {
        msg <- paste0(msg, "  + Unrooted and without node spans\n")
      } else {
        msg <- paste0(msg, "  + Unrooted, with node spans\n")
        msg <- paste0(msg, "  + PD ", signif(object@pd, 3), "\n")
      }
    } else {
      if (object@wspn) {
        msg <- paste0(msg, "  + PD ", signif(object@pd, 3), "\n")
      } else {
        msg <- paste0(msg, "  + Without node spans\n")
      }
      msg <- paste0(msg, '  + Root node is \"', object@root, '\"\n')
    }
    if (length(object@othr_slt_nms) > 0) {
      msg <- paste0(msg, "  + With additional node slots:\n")
      for (slt_nm in object@othr_slt_nms) {
        msg <- paste0(msg, "    [", slt_nm, "]\n")
      }
    }
    cat(msg)
  }
)

#' @name TreeMen-class
#' @title TreeMen-class
#' @aliases TreeMen-method
#' @param x \code{TreeMen} object
#' @param i tree index (integer or character)
#' @param object \code{TreeMen} object
#' @param max.level \code{str()} maximum level
#' @param ... additional tree objects
#' @param j missing
#' @param drop missing
#' @description S4 class for multiple phylogenetic trees
#' @slot treelst list of \code{TreeMan} objects
#' @slot ntips sum of tips per tree
#' @slot ntrees total number of trees
#' @exportClass TreeMen
#' @seealso
#' \code{\link{cTrees}}
setClass("TreeMen",
  representation = representation(
    treelst = "list", # list of TreeMan objects
    ntips = "numeric", # sum of tips per tree
    ntrees = "numeric"
  ), # number of trees in object
  validity = checkTreeMen
)

# concatenate methods
.cMenToMen <- function(treemen_1, treemen_2) {
  treelst <- c(treemen_1@treelst, treemen_2@treelst)
  treemen_1@treelst <- treelst
  treemen_1@ntips <- treemen_1@ntips + treemen_2@ntips
  treemen_1@ntrees <- treemen_1@ntrees + treemen_2@ntrees
  treemen_1
}

.cMenToMan <- function(treemen, treeman) {
  treelst <- c(treemen@treelst, treeman)
  treemen@treelst <- treelst
  treemen@ntips <- treeman@ntips + treemen@ntips
  treemen@ntrees <- treemen@ntrees + 1
  treemen
}

.cMenToAny <- function(treemen, treeobj) {
  if (class(treeobj)[1] == "TreeMan") {
    treemen <- .cMenToMan(treemen, treeobj)
  } else if (class(treeobj)[1] == "TreeMen") {
    treemen <- .cMenToMen(treemen, treeobj)
  }
  treemen
}

.cTreeObjs <- function(treemen, treeobj, ...) {
  if (nargs() > 2) {
    treemen <- .cMenToAny(treemen, treeobj)
    treemen <- .cTreeObjs(treemen, ...)
  } else {
    treemen <- .cMenToAny(treemen, treeobj)
  }
  treemen
}

#' @title cTrees
#' @description Return \code{TreeMen} of concatenated trees.
#' @details Concatenate trees into single \code{TreeMen} object.
#' @param x \code{TreeMan} or \code{TreeMen} objects
#' @param ... more \code{TreeMan} or \code{TreeMen} objects
#' @seealso
#' \code{\link{TreeMen-class}}, \code{\link{TreeMan-class}}, \code{\link{list-to-TreeMen}}
#' @examples
#'
#' trees <- cTrees(randTree(10), randTree(10))
#' @export
setGeneric("cTrees",
  signature = c("x"),
  function(x, ...) {
    standardGeneric("cTrees")
  }
)
#' @rdname TreeMan-class
#' @exportMethod cTrees
setMethod(
  "cTrees", c("TreeMan"),
  function(x, ...) {
    x <- list(x)
    x <- as(x, "TreeMen")
    x <- .cTreeObjs(x, ...)
    x
  }
)
#' @rdname TreeMen-class
#' @exportMethod cTrees
setMethod(
  "cTrees", c("TreeMen"),
  function(x, ...) {
    x <- .cTreeObjs(x, ...)
    x
  }
)

# Accessor methods
#' @rdname TreeMen-class
#' @exportMethod [[
setMethod(
  "[[", c("TreeMen", "ANY"),
  function(x, i) {
    if (!i %in% 1:length(x@treelst)) {
      stop(paste0(i, " not in tree"))
    }
    x@treelst[[i]]
  }
)
#' @rdname TreeMen-class
#' @aliases TreeMen-method
#' Extract slots from a list of trees
#' @exportMethod [
setMethod(
  "[", c("TreeMen", "character", "missing", "missing"),
  function(x, i, j, ..., drop = TRUE) {
    if (!i %in% slotNames(x)) {
      stop(paste0(i, "  not in tree"))
    }
    slot(x, i)
  }
)


# display methods
#' @rdname TreeMen-class
#' @exportMethod as.character
setMethod(
  "as.character", c("x" = "TreeMen"),
  function(x) {
    paste0("TreeMen Object of [", x@ntrees, "] trees")
  }
)
#' @rdname TreeMen-class
#' @exportMethod show
setMethod(
  "show", "TreeMen",
  function(object) {
    msg <- as.character(object)
    cat(msg)
  }
)
#' @rdname TreeMen-class
#' @exportMethod str
setMethod(
  "str", c("object" = "TreeMen"),
  function(object, max.level = 2L, ...) {
    if (is.na(max.level)) {
      stop("max.level must be numeric")
    }
    str@default(object, max.level = max.level, ...)
  }
)
#' @rdname TreeMen-class
#' @exportMethod print
setMethod(
  "print", c("x" = "TreeMen"),
  function(x) {
    msg <- as.character(x)
    print(msg)
  }
)
#' @rdname TreeMen-class
#' @exportMethod summary
setMethod(
  "summary", c("object" = "TreeMen"),
  function(object) {
    msg <- "Trees (TreeMen Object):\n"
    msg <- paste0(msg, "  + ", object@ntrees, " trees\n")
    msg <- paste0(msg, "  + ", object@ntips, " tips\n")
    cat(msg)
  }
)


#' @name pstMnp
#' @title Update prinds and tinds
#' @description Return tree with updated slots.
#' @details This function is automatically run. Only run, if you
#' are creating yor own functions to add and remove elements of the
#' \code{ndlst}.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{updateSlts}}, \code{\link{addNdmtrx}},
#' \code{\link{getAge}}
#' @export
pstMnp <- function(tree) {
  # after any adding or removing of tips and nodes,
  # these slots MUST be updated to ensure full functionality
  tree@tinds <- .getTinds(tree@ndlst)
  tree@prinds <- .getPrinds(tree@ndlst)
  tree
}

#' @name updateSlts
#' @title Update tree slots after manipulation
#' @description Return tree with updated slots.
#' @details Tree slots in the \code{TreeMan} object are usually automatically updated.
#' For certain single node manipulations they are not. Run this
#' function to update the slots.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{addNdmtrx}}, \code{\link{getAge}}
#' @export
updateSlts <- function(tree) {
  # Update the slots for a tree
  wo_pstndes <- vapply(
    tree@ndlst,
    function(n) length(n[["ptid"]]) == 0,
    logical(1)
  )
  tree@tips <- sort(names(wo_pstndes)[wo_pstndes])
  tree@ntips <- length(tree@tips)
  tree@nds <- sort(names(wo_pstndes)[!wo_pstndes])
  tree@nnds <- length(tree@nds)
  tree@all <- names(tree@ndlst)
  tree@nall <- length(tree@all)
  tree@wtxnyms <- any(vapply(
    tree@ndlst, function(n) !is.null(n[["txnym"]]),
    logical(1)
  ))
  spns <- vapply(tree@ndlst, function(n) n[["spn"]], numeric(1))
  tree@wspn <- any(spns > 0)
  if (tree@wspn) {
    tree@pd <- sum(vapply(tree@ndlst, function(n) n[["spn"]], numeric(1)))
  } else {
    tree@pd <- numeric()
  }
  tree@ply <- any(vapply(
    tree@ndlst, function(n) length(n[["ptid"]]) > 2,
    logical(1)
  ))
  tree@updtd <- TRUE
  initialize(tree)
}

#' @name addNdmtrx
#' @title Add node matrix to a tree
#' @description Return tree with node matrix added.
#' @details The node matrix makes 'enquiry'-type computations faster:
#' determining node ages, number of descendants etc. But it takes up
#' large amounts of memory and has no impact on adding or removing tips.
#' Note, trees with the node matrix can not be written to disk using the
#' 'serialization format' i.e. with \code{save} or \code{saveRDS}.
#' The matrix is generated with bigmemory's `as.big.matrix()`.
#' @param tree \code{TreeMan} object
#' @param shared T/F, should the bigmatrix be shared? See bigmemory documentation.
#' @param ... \code{as.big.matrix()} additional arguments
#' @seealso
#' \code{\link{updateSlts}}, \code{\link{rmNdmtrx}},
#' \url{https://cran.r-project.org/package=bigmemory}
#' @export
#' @examples
#' #
#' tree <- randTree(10, wndmtrx = FALSE)
#' summary(tree)
#' tree <- addNdmtrx(tree)
#' summary(tree)
addNdmtrx <- function(tree, shared = FALSE, ...) {
  if (tree@ntips < 3) {
    stop("Too small for node matrix.")
  }
  if (!checkNdlst(tree@ndlst, tree@root)) {
    stop("Invalid tree")
  }
  if (is.null(tree@ndmtrx)) {
    # generate ndmtrx
    tree@ndmtrx <- .getNdmtrxFrmLst(tree@ndlst, shared = shared, ...)
  }
  tree
}

#' @name rmNdmtrx
#' @title Remove node matrix
#' @description Return tree with memory heavy node matrix removed.
#' @details Potential uses: reduce memory load of a tree,
#' save tree using serialization methods.
#' @param tree \code{TreeMan} object
#' @seealso
#' \code{\link{addNdmtrx}}
#' @export
#' @examples
#' #
#' tree <- randTree(10)
#' summary(tree)
#' tree <- rmNdmtrx(tree)
#' summary(tree)
rmNdmtrx <- function(tree) {
  tree@ndmtrx <- NULL
  tree
}

#' @name TreeMan-to-phylo
#' @title Convert TreeMan to phylo
#' @description Return ape's \code{phylo} from a \code{TreeMan}. This will
#' preserve node labels if they are different from the default labels (n#).
#' @seealso
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples
#'
#' library(ape)
#' tree <- randTree(10)
#' tree <- as(tree, "phylo")
setAs(from = "phylo", to = "TreeMan", def = function(from, to) {
  temp_file <- tempfile(pattern = "temp_tree", fileext = ".tre")
  ape::write.tree(from, file = temp_file)
  tree <- readTree(file = temp_file)
  file.remove(temp_file)
  return(tree)
})

#' @name phylo-to-TreeMan
#' @title Convert phylo to TreeMan
#' @description Return a \code{TreeMan} from ape's \code{phylo}. This will
#' preserve node labels, if they are a alphanumeric.
#' @seealso
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples
#'
#' library(ape)
#' tree <- compute.brlen(rtree(10))
#' tree <- as(tree, "TreeMan")
setAs(from = "TreeMan", to = "phylo", def = function(from, to) {
  temp_file <- tempfile(pattern = "temp_tree", fileext = ".tre")

  # use a node label function function that writes the node labels during the
  # conversion, if they are non default.
  if (all(from["nds"] == paste0("n", seq_along(from["nds"])))) {
    nodeLabeller <- function(nd) {
      return(NULL)
    } # no node labels
  } else {
    nodeLabeller <- function(nd) {
      nd[["id"]]
    } # writes node labels
  }

  writeTree(from, file = temp_file, ndLabels = nodeLabeller)
  tree <- ape::read.tree(file = temp_file)
  file.remove(temp_file)
  return(tree)
})

#' @name multiPhylo-to-TreeMen
#' @title Convert multiPhylo to TreeMen
#' @description Return a \code{TreeMen} from ape's \code{mutlPhylo}
#' @seealso
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{TreeMen-to-multiPhylo}}
#' \code{\link{TreeMan-class}}
#' @examples
#'
#' library(ape)
#' trees <- c(rtree(10), rtree(10), rtree(10))
#' trees <- as(trees, "TreeMen")
setAs(from = "multiPhylo", to = "TreeMen", def = function(from, to) {
  temp_file <- tempfile(pattern = "temp_tree", fileext = ".tre")
  ape::write.tree(from, file = temp_file)
  tree <- readTree(file = temp_file)
  file.remove(temp_file)
  return(tree)
})

#' @name TreeMen-to-multiPhylo
#' @title Convert TreeMen to multiPhylo
#' @description Return ape's \code{multiPhylo} from a \code{TreeMen}
#' @seealso
#' \code{\link{TreeMan-to-phylo}},
#' \code{\link{phylo-to-TreeMan}},
#' \code{\link{multiPhylo-to-TreeMen}}
#' \code{\link{TreeMan-class}}
#' @examples
#'
#' library(ape)
#' trees <- cTrees(randTree(10), randTree(10), randTree(10))
#' trees <- as(trees, "multiPhylo")
setAs(from = "TreeMen", to = "multiPhylo", def = function(from, to) {
  temp_file <- tempfile(pattern = "temp_tree", fileext = ".tre")
  writeTree(from, file = temp_file)
  tree <- ape::read.tree(file = temp_file)
  file.remove(temp_file)
  return(tree)
})

#' @name list-to-TreeMen
#' @title Convert list to a TreeMen
#' @description Return a \code{TreeMen} object from a list of \code{TreeMans}
#' @seealso
#' \code{\link{TreeMen-class}}
#' @examples
#'
#' trees <- list("tree_1" = randTree(10), "tree_2" = randTree(10))
#' trees <- as(trees, "TreeMen")
# Conversion method
setAs(from = "list", to = "TreeMen", def = function(from, to) {
  ntips <- sum(unlist(lapply(from, function(tree) tree@ntips)))
  ntrees <- length(from)
  new(to, treelst = from, ntips = ntips, ntrees = ntrees)
})
