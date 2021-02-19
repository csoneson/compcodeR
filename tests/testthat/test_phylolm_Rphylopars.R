context("Rphylopars vs phylolm")

test_that("Rphylopars vs phylolm", {
  skip_if_not_installed("Rphylopars")
  skip_if_not_installed("phylolm")
  
  ## Tree
  set.seed(1289)
  ntips <- 20
  tree <- ape::rphylo(ntips, 0.1, 0)
  
  ## Replicates
  r <- 3
  tree_rep <- add_replicates(tree, r)
  
  ## traits
  sigma2_phylo <- 1
  sigma2_intra <- 0.1
  resids <- phylolm::rTrait(n = 1,
                            phy = tree_rep,
                            model = "BM",
                            parameters = list(sigma2 = sigma2_phylo),
                            plot.tree = FALSE)
  resids <- resids + rnorm(ntips * r, mean = 0, sd = sqrt(sigma2_intra))
  traits <- data.frame(trait = resids)
  
  ##############################################################################
  ### BM
  
  ## phylolm
  fit_phylolm <- phylolm::phylolm(trait~1, traits, tree_rep, measurement_error = TRUE)
  fit_phylolm_lambda <- phylolm::phylolm(trait~1, traits, tree_rep, model = "lambda")
  
  ## Test that lambda is the same as measurement error
  lambda_error <- fit_phylolm$sigma2 / (fit_phylolm$sigma2_error / max(ape::vcv(tree_rep)) + fit_phylolm$sigma2)
  expect_equal(lambda_error, fit_phylolm_lambda$optpar, tolerance = 1e-5)
  
  ## Rphylopars
  traits$species <- sub("\\_.", "", rownames(traits))
  traits <- traits[, c(2, 1)]
  fit_phylopars <- Rphylopars::phylopars(traits, tree, REML = FALSE)
  
  ## Test phylo variance
  expect_equivalent(fit_phylopars$pars$phylocov, fit_phylolm$sigma2, tolerance = 1e-6)
  ## Test pheno variance
  expect_equivalent(fit_phylopars$pars$phenocov, fit_phylolm$sigma2_error, tolerance = 1e-6)
  ## Likelihood
  expect_equivalent(fit_phylopars$logLik, fit_phylolm$logLik)
  
  ##############################################################################
  ### OU
  
  ## Rphylopars
  fit_phylopars <- Rphylopars::phylopars(traits, tree, model = "OU", REML = FALSE)
  
  ## fiddle with the tree so that there is no exactly zero branch lengths.
  tree_rep_per <- tree_rep
  tree_rep_per$edge.length <- tree_rep_per$edge.length + min(tree_rep$edge.length[tree_rep$edge.length != 0]) / 100
  tree_rep_per <- phangorn::nnls.tree(ape::cophenetic.phylo(tree_rep_per), tree_rep_per, rooted = TRUE, trace = 0)
  
  ## phylolm
  fit_phylolm <- phylolm::phylolm(trait ~ 1, traits, tree_rep_per,
                                  model = "OUfixedRoot", measurement_error = TRUE)
  
  ## Test phylo variance
  expect_equivalent(fit_phylopars$pars$phylocov, fit_phylolm$sigma2, tolerance = 1e-6)
  ## Test pheno variance
  expect_equivalent(fit_phylopars$pars$phenocov, fit_phylolm$sigma2_error, tolerance = 1e-6)
  ## Test alpha
  expect_equivalent(fit_phylopars$model$alpha, fit_phylolm$optpar, tolerance = 1e-6)
  ## Likelihood
  expect_equivalent(fit_phylopars$logLik, fit_phylolm$logLik)
  
  ##############################################################################
  ### OU - lambda transform
  
  ## OU
  fit_ou <- phylolm::phylolm(trait ~ 1, traits, tree_rep_per, model = "OUfixedRoot", measurement_error = TRUE)
  
  ## Transform the tree
  tree_ou <- phylolm::transf.branch.lengths(tree_rep_per, "OUfixedRoot", parameters = list(alpha = fit_ou$optpar))
  tree_ou <- tree_ou$tree
  
  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_ou, model = "lambda", measurement_error = FALSE)
  
  ## lambda on OU
  tilde_t <- max(ape::vcv(tree_ou)) / (2 * fit_ou$optpar)
  lambda_ou_error <- fit_ou$sigma2 * tilde_t / (fit_ou$sigma2_error + fit_ou$sigma2 * tilde_t)
  
  ## Both are equal
  expect_equivalent(lambda_ou_error, fit_lambda$optpar, tolerance = 1e-5)
  
  ## Variances
  expect_equivalent(fit_ou$sigma2 / (2 * fit_ou$optpar) + fit_ou$sigma2_error / max(ape::vcv(tree_ou)), fit_lambda$sigma2, tolerance = 1e-4)
  
  ## Transform again
  tree_ou_lambda <- phylolm::transf.branch.lengths(tree_ou, "lambda", parameters = list(lambda = lambda_ou_error))
  tree_ou_lambda <- tree_ou_lambda$tree
  
  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_ou_lambda, model = "BM", measurement_error = FALSE)
  
  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_ou$logLik)
  
  ## Variances
  expect_equivalent(fit_bm$sigma2, fit_lambda$sigma2, tolerance = 1e-4)
  
  ##############################################################################
  ### delta - lambda transform
  
  ## delta
  fit_delta <- phylolm::phylolm(trait ~ 1, traits, tree_rep_per, model = "delta", measurement_error = TRUE)
  
  ## Transform the tree
  tree_delta <- phylolm::transf.branch.lengths(tree_rep_per, "delta", parameters = list(delta = fit_delta$optpar))
  tree_delta <- tree_delta$tree
  
  ## Lambda on the transformed tree
  fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_delta, model = "lambda", measurement_error = FALSE)
  
  ## lambda on delta
  tilde_t <- max(ape::vcv(tree_delta))
  lambda_delta_error <- fit_delta$sigma2 * tilde_t / (fit_delta$sigma2_error + fit_delta$sigma2 * tilde_t)
  
  ## Both are equal
  expect_equivalent(lambda_delta_error, fit_lambda$optpar, tolerance = 1e-4)
  
  ## Transform again
  tree_delta_lambda <- phylolm::transf.branch.lengths(tree_delta, "lambda", parameters = list(lambda = lambda_delta_error))
  tree_delta_lambda <- tree_delta_lambda$tree
  
  ## Fit BM on transformed tree
  fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_delta_lambda, model = "BM", measurement_error = FALSE)
  
  ## Same likelihood
  expect_equal(fit_bm$logLik, fit_delta$logLik)
  
  # ##############################################################################
  # ### EB - lambda transform
  # 
  # ## EB
  # fit_EB <- phylolm::phylolm(trait ~ 1, traits, tree_rep_per, model = "EB", measurement_error = TRUE)
  # 
  # ## Transform the tree
  # tree_EB <- phylolm::transf.branch.lengths(tree_rep_per, "EB", parameters = list(EB = fit_EB$optpar))
  # tree_EB <- tree_EB$tree
  # 
  # ## Lambda on the transformed tree
  # fit_lambda <- phylolm::phylolm(trait ~ 1, traits, tree_EB, model = "lambda", measurement_error = FALSE)
  # 
  # ## lambda on EB
  # tilde_t <- max(ape::vcv(tree_EB))
  # lambda_EB_error <- fit_EB$sigma2 * tilde_t / (fit_EB$sigma2_error + fit_EB$sigma2 * tilde_t)
  # 
  # ## Both are equal
  # expect_equivalent(lambda_EB_error, fit_lambda$optpar, tolerance = 1e-4)
  # 
  # ## Transform again
  # tree_EB_lambda <- phylolm::transf.branch.lengths(tree_EB, "lambda", parameters = list(lambda = lambda_EB_error))
  # tree_EB_lambda <- tree_EB_lambda$tree
  # 
  # ## Fit BM on transformed tree
  # fit_bm <- phylolm::phylolm(trait ~ 1, traits, tree_EB_lambda, model = "BM", measurement_error = FALSE)
  # 
  # ## Same likelihood
  # expect_equal(fit_bm$logLik, fit_EB$logLik)
  
})

test_that("phylolm p-values", {
  skip_if_not_installed("Rphylopars")
  skip_if_not_installed("phylolm")
  
  ## Tree
  set.seed(1289)
  ntips <- 10
  tree <- ape::rphylo(ntips, 0.1, 0)
  
  ## Replicates
  r <- 2
  tree_rep <- add_replicates(tree, r)
  
  ## traits
  sigma2_phylo <- 1
  sigma2_intra <- 0.1
  resids <- phylolm::rTrait(n = 1,
                            phy = tree_rep,
                            model = "BM",
                            parameters = list(sigma2 = sigma2_phylo),
                            plot.tree = FALSE)
  resids <- resids + rnorm(ntips * r, mean = 0, sd = sqrt(sigma2_intra))
  x <- sample(c(0, 1), ntips * r, replace = TRUE)
  traits <- data.frame(trait = resids + 0.01 * x, x = x)
  
  ## phylolm
  fit_phylolm <- phylolm::phylolm(trait ~ x, traits, tree_rep, measurement_error = FALSE)
  fit_phylolm_error <- phylolm::phylolm(trait ~ x, traits, tree_rep, measurement_error = TRUE)
  fit_phylolm_lambda <- phylolm::phylolm(trait ~ x, traits, tree_rep, model = "lambda")
  
  summary(fit_phylolm)
  expect_equal(summary(fit_phylolm_error)$coefficients[, "p.value"],
               summary(fit_phylolm_lambda)$coefficients[, "p.value"],
               tolerance = 1e-3)
})



