context("Lengths and length factors")

test_that("computeFactorLengths is equivalent to trueProba", {
  set.seed(20200430)
  
  n.sample <- 3
  n.vars <- 5
  
  ## Mock data
  count.matrix <- matrix(c(5,3,8,17,23,42,10,13,27,752,615,1203,1507,1225,2455), ncol = n.sample, byrow = T)
  truemeans <- rowMeans(count.matrix)
  leng <- matrix(c(rep(1000, 2 * n.vars - 1), 2000, 1:n.vars*1000), ncol = n.sample, byrow = F)
  
  ## True proba
  trueMeansLeng <- diag(truemeans) %*% leng
  trueProbas <- trueMeansLeng %*% diag(1 / colSums(trueMeansLeng))
  
  ## nfact
  lfact <- computeFactorLengths(leng, truemeans, sum(truemeans))
  
  ## Equal
  expect_equal(trueProbas, diag(truemeans / sum(truemeans)) %*% lfact)
  
  ## Sequencing depths
  seqdepth <- 2500
  nfacts <- runif(n.sample, min = 0.7, max = 1.4)
  seq.depths <- nfacts * seqdepth
  
  ## expectations
  exp_1 <- trueProbas %*% diag(seq.depths)
  expect_equal(colSums(exp_1), seq.depths)
  
  exp_2 <- diag(truemeans / sum(truemeans)) %*% lfact %*% diag(seq.depths)
  expect_equal(exp_1, exp_2)
})

test_that("computeFactorLengths is equivalent to trueProba", {
  set.seed(20200430)
  
  ntaxa <- 10
  nrep <- 3
  n.sample <- ntaxa * nrep
  n.vars <- 50
  
  ## Mock tree
  set.seed("17900714")
  tree <- ape::rphylo(ntaxa, birth = 0.1, death = 0)
  # rescale to unit height
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  # rename tims
  tree$tip.label <- paste0("t", 1:ntaxa)
  tree <- add_replicates(tree, nrep)
  
  id.species <- rep(1:ntaxa, each = nrep)
  names(id.species) <- tree$tip.label
  id.species <- as.factor(id.species)
  
  id.species <- checkSpecies(id.species, "id.species", tree, 1e-10, TRUE)
  
  ## Mock data
  count.matrix <- matrix(rpois(n.sample * n.vars, 1:10 * 10), ncol = n.sample, byrow = T)
  truemeans <- rowMeans(count.matrix)
  leng <- matrix(c(rep(1000, 2 * n.vars - 1), 2000, 1:n.vars*1000), ncol = n.sample, byrow = F)
  lengths.relmeans <- 1:n.vars * 1000
  lengths.dispersions <- 1:n.vars / 10
  leng <- generateLengths(id.species, lengths.relmeans, lengths.dispersions)
  
  expect_equal(leng[, 1], leng[, 2])
  expect_equal(leng[, 4], leng[, 5])
  
  ## True proba
  trueMeansLeng <- diag(truemeans) %*% leng
  trueProbas <- trueMeansLeng %*% diag(1 / colSums(trueMeansLeng))
  
  ## nfact
  lfact <- computeFactorLengths(leng, truemeans, sum(truemeans))
  
  ## Equal
  expect_equal(trueProbas, diag(truemeans / sum(truemeans)) %*% lfact)
  
  ## Sequencing depths
  seqdepth <- 2500
  nfacts <- runif(n.sample, min = 0.7, max = 1.4)
  seq.depths <- nfacts * seqdepth
  
  ## expectations
  exp_1 <- trueProbas %*% diag(seq.depths)
  expect_equal(colSums(exp_1), seq.depths)
  
  exp_2 <- diag(truemeans / sum(truemeans)) %*% lfact %*% diag(seq.depths)
  expect_equal(exp_1, exp_2)
})