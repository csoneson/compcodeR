context("Negative Binomial to Poisson Log Normal")

test_that("NB to PLN simple", {
  set.seed(18420318)
  
  n <- 20000000
  
  ## NB
  mean_nb <- 500
  dispersion_nb <- 0.2
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- rnbinom(n = n,
                       mu = mean_nb, 
                       size = 1 / dispersion_nb)
  
  ## PLN
  params_PLN <- NB_to_PLN(mean_nb, dispersion_nb)
  
  sample_ln <- rnorm(n = n,
                      mean = params_PLN$log_means_pln,
                      sd = sqrt(params_PLN$log_variances_pln))
  sample_pln <- rpois(n, exp(sample_ln))
  
  ## Comparison
  expect_equal(mean(sample_nb) - mean_nb, 0, tolerance = 0.05)
  expect_equal(sd(sample_nb) - sd_nb, 0, tolerance = 0.05)
  
  expect_equal(mean(sample_pln) - mean_nb, 0, tolerance = 0.05)
  expect_equal(sd(sample_pln) - sd_nb, 0, tolerance = 0.07)
})

test_that("NB to PLN phylo - star tree", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phylolm")
  
  set.seed(18420318)
  
  ## Parameters
  n <- 200000
  ntaxa <- 4
  
  ## Tree
  tree <- ape::stree(ntaxa, type = "star")
  tree <- ape::compute.brlen(tree, runif, min = 0, max = 1)
  tree <- phytools::force.ultrametric(tree) 
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  
  ## NB
  mean_nb <- 1:ntaxa * 100
  dispersion_nb <- 1:ntaxa/2 / 100
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- t(matrix(rnbinom(n = ntaxa * n,
                                mu = mean_nb, 
                                size = 1 / dispersion_nb), nrow = ntaxa))
  
  ## PLN
  names(mean_nb) <- tree$tip.label
  names(dispersion_nb) <- tree$tip.label
  
  params_PLN <- get_poisson_log_normal_parameters(rep(1, n) %*% t(mean_nb), rep(1, n) %*% t(dispersion_nb), 1.0)
  
  sample_ppln <- simulatePhyloPoissonLogNormal(tree, params_PLN$log_means, params_PLN$log_variance_phylo, params_PLN$log_variance_sample)
  
  sample_ln <- sample_ppln$log_lambda
  sample_pln <- sample_ppln$counts
  rm(sample_ppln)
  
  mean_ln <- params_PLN$log_means[1, ]
  sd_ln <- sqrt((params_PLN$log_variance_phylo + params_PLN$log_variance_sample)[1, ])
  
  ## Comparisons NB
  expect_equivalent(colMeans(sample_nb) - mean_nb, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_nb) - sd_nb, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparison log lambda
  expect_equivalent(colMeans(sample_ln) - mean_ln, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_ln) - sd_ln, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparisons PLN
  expect_equivalent(colMeans(sample_pln) - mean_nb, rep(0, ntaxa), tolerance = 0.1)
  expect_equivalent(matrixStats::colSds(sample_pln) - sd_nb, rep(0, ntaxa), tolerance = 0.1)
})

test_that("NB to PLN phylo - random tree", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phylolm")
  
  set.seed(18420318)
  
  ## Parameters
  n <- 1000000
  ntaxa <- 4
  
  ## Tree
  tree <- ape::rtree(ntaxa)
  tree <- ape::compute.brlen(tree, runif, min = 0, max = 1)
  tree <- phytools::force.ultrametric(tree) 
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  
  ## NB
  mean_nb <- 1:ntaxa * 100
  dispersion_nb <- 1:ntaxa/2 / 100
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- t(matrix(rnbinom(n = ntaxa * n,
                                mu = mean_nb, 
                                size = 1 / dispersion_nb), nrow = ntaxa))
  
  ## PLN
  names(mean_nb) <- tree$tip.label
  names(dispersion_nb) <- tree$tip.label
  
  params_PLN <- get_poisson_log_normal_parameters(rep(1, n) %*% t(mean_nb), rep(1, n) %*% t(dispersion_nb), 1.0)
  
  sample_ppln <- simulatePhyloPoissonLogNormal(tree, params_PLN$log_means, params_PLN$log_variance_phylo, params_PLN$log_variance_sample)
  
  sample_ln <- sample_ppln$log_lambda
  sample_pln <- sample_ppln$counts
  rm(sample_ppln)
  
  mean_ln <- params_PLN$log_means[1, ]
  sd_ln <- sqrt((params_PLN$log_variance_phylo + params_PLN$log_variance_sample)[1, ])
  
  ## Comparisons NB
  expect_equivalent(colMeans(sample_nb) - mean_nb, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_nb) - sd_nb, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparison log lambda
  expect_equivalent(colMeans(sample_ln) - mean_ln, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_ln) - sd_ln, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparisons PLN
  expect_equivalent(colMeans(sample_pln) - mean_nb, rep(0, ntaxa), tolerance = 0.1)
  expect_equivalent(matrixStats::colSds(sample_pln) - sd_nb, rep(0, ntaxa), tolerance = 0.1)
  
  ## Phylogenetic covariances
  C_tree <- ape::vcv(tree)
  V_ln <- C_tree * params_PLN$log_variance_phylo[1] + diag(params_PLN$log_variance_sample[1, ])
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_ln[i, j] != 0) {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]), V_ln[i, j], tolerance = 0.0001, scale = V_ln[i, j])
      } else {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]) - V_ln[i, j], 0, tolerance = 0.0001)
      }
    }
  }
  
  ## Phylogenetic covariances for Counts
  V_tot <- (exp(V_ln) - 1) * mean_nb %*% t(mean_nb)
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_tot[i, j] != 0) {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]), V_tot[i, j], tolerance = 0.1, scale = V_tot[i, j])
      } else {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]) - V_tot[i, j], 0, tolerance = 0.4)
      }
    }
  }
})

test_that("NB to PLN phylo - star tree with repetitions", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phylolm")
  
  set.seed(18420318)
  
  ## Parameters
  n <- 3000000
  ntaxa <- 2
  
  ## Tree
  tree <- ape::stree(ntaxa, type = "star")
  tree <- ape::compute.brlen(tree, runif, min = 0, max = 1)
  tree <- phytools::force.ultrametric(tree) 
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  
  ## Repetitions
  r <- 2
  ntaxa <- r * ntaxa
  
  tree_rep <- tree
  # Add replicates
  for (tip_label in tree$tip.label) {
    for (rep in 1:r) {
      tree_rep <- phytools::bind.tip(tree_rep, tip.label = paste0(tip_label, "_", rep),
                                     where = which(tree_rep$tip.label == tip_label))
    }
  }
  # Remove original tips
  tree_rep <- ape::drop.tip(tree_rep, tree$tip.label)
  tree <- tree_rep
  
  ## NB
  mean_nb <- 1:ntaxa * 100
  dispersion_nb <- 1:ntaxa/2 / 100
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- t(matrix(rnbinom(n = ntaxa * n,
                                mu = mean_nb, 
                                size = 1 / dispersion_nb), nrow = ntaxa))
  
  ## PLN
  names(mean_nb) <- tree$tip.label
  names(dispersion_nb) <- tree$tip.label
  
  params_PLN <- get_poisson_log_normal_parameters(rep(1, n) %*% t(mean_nb), rep(1, n) %*% t(dispersion_nb), 1.0)
  
  sample_ppln <- simulatePhyloPoissonLogNormal(tree, params_PLN$log_means, params_PLN$log_variance_phylo, params_PLN$log_variance_sample)
  
  sample_ln <- sample_ppln$log_lambda
  sample_pln <- sample_ppln$counts
  rm(sample_ppln)
  
  mean_ln <- params_PLN$log_means[1, ]
  sd_ln <- sqrt((params_PLN$log_variance_phylo + params_PLN$log_variance_sample)[1, ])
  
  ## Comparisons NB
  expect_equivalent(colMeans(sample_nb) - mean_nb, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_nb) - sd_nb, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparison log lambda
  expect_equivalent(colMeans(sample_ln) - mean_ln, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_ln) - sd_ln, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparisons PLN
  expect_equivalent(colMeans(sample_pln) - mean_nb, rep(0, ntaxa), tolerance = 0.1)
  expect_equivalent(matrixStats::colSds(sample_pln) - sd_nb, rep(0, ntaxa), tolerance = 0.1)
  
  ## Phylogenetic covariances
  C_tree <- ape::vcv(tree)
  V_ln <- C_tree * params_PLN$log_variance_phylo[1] + diag(params_PLN$log_variance_sample[1, ])
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_ln[i, j] != 0) {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]), V_ln[i, j], tolerance = 0.0001, scale = V_ln[i, j])
      } else {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]) - V_ln[i, j], 0, tolerance = 0.0001)
      }
    }
  }
  
  ## Phylogenetic covariances for Counts
  V_tot <- (exp(V_ln) - 1) * mean_nb %*% t(mean_nb)
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_tot[i, j] != 0) {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]), V_tot[i, j], tolerance = 0.05, scale = V_tot[i, j])
      } else {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]) - V_tot[i, j], 0, tolerance = 0.6)
      }
    }
  }
})

test_that("NB to PLN phylo - random tree - OU", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phylolm")
  
  set.seed(18420318)
  
  ## Parameters
  n <- 3000000
  ntaxa <- 4
  selection.strength <- 1
  
  ## Tree
  tree <- ape::rtree(ntaxa)
  tree <- ape::compute.brlen(tree, runif, min = 0, max = 1)
  tree <- phytools::force.ultrametric(tree) 
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  
  ## NB
  mean_nb <- 1:ntaxa * 100
  dispersion_nb <- 1:ntaxa/2 / 100
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- t(matrix(rnbinom(n = ntaxa * n,
                                mu = mean_nb, 
                                size = 1 / dispersion_nb), nrow = ntaxa))
  
  ## PLN
  names(mean_nb) <- tree$tip.label
  names(dispersion_nb) <- tree$tip.label
  
  params_PLN <- get_poisson_log_normal_parameters(rep(1, n) %*% t(mean_nb), rep(1, n) %*% t(dispersion_nb), 1.0)
  
  sample_ppln <- simulatePhyloPoissonLogNormal(tree,
                                               params_PLN$log_means,
                                               params_PLN$log_variance_phylo,
                                               params_PLN$log_variance_sample,
                                               model_process = "OU",
                                               selection.strength = selection.strength)
  
  sample_ln <- sample_ppln$log_lambda
  sample_pln <- sample_ppln$counts
  rm(sample_ppln)
  
  mean_ln <- params_PLN$log_means[1, ]
  sd_ln <- sqrt((params_PLN$log_variance_phylo + params_PLN$log_variance_sample)[1, ])
  
  ## Comparisons NB
  expect_equivalent(colMeans(sample_nb) - mean_nb, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_nb) - sd_nb, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparison log lambda
  expect_equivalent(colMeans(sample_ln) - mean_ln, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_ln) - sd_ln, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparisons PLN
  expect_equivalent(colMeans(sample_pln) - mean_nb, rep(0, ntaxa), tolerance = 0.1)
  expect_equivalent(matrixStats::colSds(sample_pln) - sd_nb, rep(0, ntaxa), tolerance = 0.1)
  
  ## Phylogenetic covariances
  C_tree <- ape::vcv(phylolm::transf.branch.lengths(tree, model = "OUfixedRoot", parameters = list(alpha = selection.strength))$tree)
  V_ln <- C_tree * params_PLN$log_variance_phylo[1] / -expm1(-2 * selection.strength) + diag(params_PLN$log_variance_sample[1, ])
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_ln[i, j] != 0) {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]), V_ln[i, j], tolerance = 0.00001, scale = V_ln[i, j])
      } else {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]) - V_ln[i, j], 0, tolerance = 0.0001)
      }
    }
  }
  
  ## Phylogenetic covariances for Counts
  V_tot <- (exp(V_ln) - 1) * mean_nb %*% t(mean_nb)
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_tot[i, j] != 0) {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]), V_tot[i, j], tolerance = 0.06, scale = V_tot[i, j])
      } else {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]) - V_tot[i, j], 0, tolerance = 0.6)
      }
    }
  }
})

test_that("NB to PLN phylo - random tree - Not Unit Length", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phylolm")
  
  set.seed(18420318)
  
  ## Parameters
  n <- 4000000
  ntaxa <- 4
  selection.strength <- 1
  
  ## Tree
  tree <- ape::rtree(ntaxa)
  tree <- ape::compute.brlen(tree, runif, min = 0, max = 1)
  tree <- phytools::force.ultrametric(tree) 
  tree_height <- ape::vcv(tree)[1, 1]

  ## NB
  mean_nb <- 1:ntaxa * 100
  dispersion_nb <- 1:ntaxa/2 / 100
  
  sd_nb <- sqrt(mean_nb + dispersion_nb * mean_nb^2)
  
  sample_nb <- t(matrix(rnbinom(n = ntaxa * n,
                                mu = mean_nb, 
                                size = 1 / dispersion_nb), nrow = ntaxa))
  
  ## PLN
  names(mean_nb) <- tree$tip.label
  names(dispersion_nb) <- tree$tip.label
  
  params_PLN <- get_poisson_log_normal_parameters(rep(1, n) %*% t(mean_nb), rep(1, n) %*% t(dispersion_nb), 1.0)
  
  sample_ppln <- simulatePhyloPoissonLogNormal(tree,
                                               params_PLN$log_means,
                                               params_PLN$log_variance_phylo,
                                               params_PLN$log_variance_sample)
  
  sample_ln <- sample_ppln$log_lambda
  sample_pln <- sample_ppln$counts
  rm(sample_ppln)
  
  mean_ln <- params_PLN$log_means[1, ]
  sd_ln <- sqrt((params_PLN$log_variance_phylo + params_PLN$log_variance_sample)[1, ])
  
  ## Comparisons NB
  expect_equivalent(colMeans(sample_nb) - mean_nb, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_nb) - sd_nb, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparison log lambda
  expect_equivalent(colMeans(sample_ln) - mean_ln, rep(0, ntaxa), tolerance = 0.05)
  expect_equivalent(matrixStats::colSds(sample_ln) - sd_ln, rep(0, ntaxa), tolerance = 0.05)
  
  ## Comparisons PLN
  expect_equivalent(colMeans(sample_pln) - mean_nb, rep(0, ntaxa), tolerance = 0.1)
  expect_equivalent(matrixStats::colSds(sample_pln) - sd_nb, rep(0, ntaxa), tolerance = 0.1)
  
  ## Phylogenetic covariances
  C_tree <- ape::vcv(tree)
  V_ln <- C_tree * params_PLN$log_variance_phylo[1] / tree_height + diag(params_PLN$log_variance_sample[1, ])
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_ln[i, j] != 0) {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]), V_ln[i, j], tolerance = 0.0001, scale = V_ln[i, j])
      } else {
        expect_equivalent(cov(sample_ln[, i], sample_ln[, j]) - V_ln[i, j], 0, tolerance = 0.0001)
      }
    }
  }
  
  ## Phylogenetic covariances for Counts
  V_tot <- (exp(V_ln) - 1) * mean_nb %*% t(mean_nb)
  for (i in 1:(ntaxa-1)) {
    for (j in (i+1):ntaxa) {
      if (V_tot[i, j] != 0) {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]), V_tot[i, j], tolerance = 0.06, scale = V_tot[i, j])
      } else {
        expect_equivalent(cov(sample_pln[, i], sample_pln[, j]) - V_tot[i, j], 0, tolerance = 0.6)
      }
    }
  }
})
