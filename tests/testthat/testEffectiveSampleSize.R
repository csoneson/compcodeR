context("Effective Sample Size")

test_that("Errors", {
  
  ##############################################################################
  ### Errors
  ntips <- 20
  id_cond_3 <- sample(0:2, ntips, replace = TRUE)
  expect_error(nEffPhylolm(NULL, id_cond_3, "BM", 0), "only two conditions")
  
})

test_that("Naive vs phylolm", {
  skip_if_not_installed("phylolm")
  
  ## Tree
  set.seed(1289)
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)
  id_cond <- sample(c(0, 1), ntips, replace = TRUE)
  names(id_cond) <- tree$tip.label
  
  ## Replicates
  r <- 3
  tree_rep <- add_replicates(tree, r)
  id_cond_rep <- rep(id_cond, each = r)
  names(id_cond_rep) <- as.vector(sapply(names(id_cond), function(x) paste(x, 1:3, sep = "_")))
  
  ##############################################################################
  ### BM
  
  expect_equal(nEffNaive(tree, id_cond, "BM", 0),
               nEffPhylolm(tree, id_cond, "BM", 0))
  
  expect_equal(nEffNaive(tree, id_cond, "BM", 0),
               nEffSchur(tree, id_cond, "BM", 0))
  
  # microbenchmark::microbenchmark(nEffNaive(tree, id_cond, "BM", 0),
  #                                nEffPhylolm(tree, id_cond, "BM", 0),
  #                                nEffSchur(tree, id_cond, "BM", 0))
  # The three are equivalent.
  
  expect_equal(tree$edge.length, prune_tree_one_obs(tree_rep)$edge.length)

  expect_warning(expect_equal(nEffNaive(tree_rep, id_cond_rep, "BM", 0),
                              nEffNaive(tree, id_cond, "BM", 0)),
                              "not sorted in the correct order")
  
  ##############################################################################
  ### OU
  
  expect_equal(nEffNaive(tree, id_cond, "OU", 3),
               nEffPhylolm(tree, id_cond, "OU", 3))
  
  expect_equal(nEffNaive(tree, id_cond, "OU", 3),
               nEffSchur(tree, id_cond, "OU", 3))
  
  # microbenchmark::microbenchmark(nEffNaive(tree, id_cond, "OU", 3),
  #                                nEffPhylolm(tree, id_cond, "OU", 3),
  #                                nEffSchur(tree, id_cond, "OU", 3))
  
  expect_equal(tree$edge.length, prune_tree_one_obs(tree_rep)$edge.length)
  
  expect_warning(expect_equal(nEffNaive(tree_rep, id_cond_rep, "OU", 3),
                              nEffNaive(tree, id_cond, "OU", 3)),
                 "not sorted in the correct order")
  
})

test_that("nEff Ratio", {
  skip_if_not_installed("phylolm")
  
  ## Trees
  set.seed(1289)
  ntips <- 2^5
  tree <- ape::stree(ntips, "balanced")
  star_tree <- ape::stree(ntips, "star")
  tree <- ape::compute.brlen(tree)
  star_tree <- ape::compute.brlen(star_tree)
  
  ## Alt cond : nEff greater than 1
  id_cond <- rep(rep(0:1, each = 2), ntips / 4)
  names(id_cond) <- tree$tip.label
  expect_true(nEffRatio(tree, id_cond, "BM", 0) > 1.0)
  expect_equal(nEffRatio(star_tree, id_cond, "BM", 0), 1.0)
  expect_equal(nEffRatio(NULL, id_cond, "BM", 0), 1.0)
  
  ## Bloc cond : nEff smaller than 1
  id_cond <- rep(0:1, each = ntips / 2)
  names(id_cond) <- tree$tip.label
  expect_true(nEffRatio(tree, id_cond, "BM", 0) < 1.0)
  expect_equal(nEffRatio(star_tree, id_cond, "BM", 0), 1.0)
  expect_equal(nEffRatio(NULL, id_cond, "BM", 0), 1.0)
  
})



