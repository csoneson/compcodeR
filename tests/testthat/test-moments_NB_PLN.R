context("Moments of the Negative Binomial and Poisson Log Normal are the same")

test_that("Moments of the Negative Binomial and Log Normal are the same", {
  
  tdir <- tempdir()
  
  ## Tree and conds
  n <- 20
  tree <- ape::stree(n)
  tree$edge.length <- rep(10, nrow(tree$edge))
  
  id_species <- tree$tip.label
  id_species <- factor(id_species)
  names(id_species) <- tree$tip.label
  
  ## Condn alternate tips
  species_names <- tree$tip.label
  cond_species <- rep(c(1, 2), length(species_names) / 2)
  names(cond_species) <- species_names
  
  id_cond <- id_species
  id_cond <- cond_species[as.vector(id_cond)]
  id_cond <- as.factor(id_cond)
  
  ## Parameters
  n.vars <- 1000
  relmean <- 100
  reldisp <- 0.1
  seqdepth <- 10000
  
  true_mean <- seqdepth / n.vars
  true_var <- reldisp * true_mean^2 + true_mean
  true_sd <- sqrt(true_var)
  
  ## No tree
  set.seed(18570823)
  dataset <- "no_tree"
  no_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = 0,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = 1,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  ## us tree
  set.seed(18570823)
  dataset <- "tree"
  us_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = 0,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = 1,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       tree = tree,
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  

  ## Comparisons
  mean_no_tree <- mean(no_tree_sim@count.matrix)
  sd_no_tree <- sd(no_tree_sim@count.matrix)
  
  mean_us_tree <- mean(us_tree_sim@count.matrix)
  sd_us_tree <- sd(us_tree_sim@count.matrix)
  
  testthat::expect_equal(mean_no_tree, true_mean, tol = 1e-2)
  testthat::expect_equal(mean_us_tree, true_mean, tol = 1e-2)
  
  testthat::expect_equal(sd_no_tree, true_sd, tol = 1e-2)
  testthat::expect_equal(sd_us_tree, true_sd, tol = 1e-2)
})

test_that("Moments of the Negative Binomial and Log Normal are the same, with shift", {
  
  tdir <- tempdir()
  
  ## Tree and conds
  n <- 20
  tree <- ape::stree(n)
  tree$edge.length <- rep(10, nrow(tree$edge))
  
  id_species <- tree$tip.label
  id_species <- factor(id_species)
  names(id_species) <- tree$tip.label
  
  ## Condn alternate tips
  species_names <- tree$tip.label
  cond_species <- rep(c(1, 2), length(species_names) / 2)
  names(cond_species) <- species_names
  
  id_cond <- id_species
  id_cond <- cond_species[as.vector(id_cond)]
  id_cond <- as.factor(id_cond)
  
  ## Parameters
  n.vars <- 5000
  relmean <- 100
  reldisp <- 0.1
  seqdepth <- 50000
  
  n.diffexp <- n.vars / 2
  effect.size <- 10
  
  true_mean_1 <- seqdepth / n.vars
  true_mean_2 <- 2 * seqdepth * (effect.size + 1) / n.vars / (effect.size + 2)
  true_var_1 <- reldisp * true_mean_1^2 + true_mean_1
  true_sd_1 <- sqrt(true_var_1)
  true_var_2 <- reldisp * true_mean_2^2 + true_mean_2
  true_sd_2 <- sqrt(true_var_2)
  
  ## No tree
  set.seed(18570823)
  dataset <- "no_tree"
  no_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  ## us tree
  set.seed(18570823)
  dataset <- "tree"
  us_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       tree = tree,
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  
  ## Comparisons
  mean_no_tree <- colMeans(no_tree_sim@count.matrix[1:n.diffexp, ])
  sd_no_tree <- matrixStats::colSds(no_tree_sim@count.matrix[1:n.diffexp, ])
  
  mean_us_tree <- colMeans(us_tree_sim@count.matrix[1:n.diffexp, ])
  sd_us_tree <-  matrixStats::colSds(us_tree_sim@count.matrix[1:n.diffexp, ])
  
  testthat::expect_equal(mean(mean_no_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_no_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(mean_us_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_us_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  testthat::expect_equal(mean(sd_us_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 2]),
                         mean(sd_us_tree[id_cond == 2]), tol = 1e-2)
  
})

test_that("Moments of the Negative Binomial and Log Normal are the same, with shift, random tree", {
  
  tdir <- tempdir()
  
  ## Tree and conds
  n <- 20
  set.seed(1789)
  tree <- tree <- ape::rphylo(n, 1, 0)

  id_species <- tree$tip.label
  id_species <- factor(id_species)
  names(id_species) <- tree$tip.label
  
  ## Condn alternate tips
  species_names <- tree$tip.label
  cond_species <- rep(c(1, 2), length(species_names) / 2)
  names(cond_species) <- species_names
  
  id_cond <- id_species
  id_cond <- cond_species[as.vector(id_cond)]
  id_cond <- as.factor(id_cond)
  
  ## Parameters
  n.vars <- 5000
  relmean <- 100
  reldisp <- 0.1
  seqdepth <- 50000
  
  n.diffexp <- n.vars / 2
  effect.size <- 10
  
  true_mean_1 <- seqdepth / n.vars
  true_mean_2 <- 2 * seqdepth * (effect.size + 1) / n.vars / (effect.size + 2)
  true_var_1 <- reldisp * true_mean_1^2 + true_mean_1
  true_sd_1 <- sqrt(true_var_1)
  true_var_2 <- reldisp * true_mean_2^2 + true_mean_2
  true_sd_2 <- sqrt(true_var_2)
  
  ## No tree
  set.seed(18570823)
  dataset <- "no_tree"
  no_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  ## us tree
  set.seed(18570823)
  dataset <- "tree"
  us_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       tree = tree,
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  
  ## Comparisons
  mean_no_tree <- colMeans(no_tree_sim@count.matrix[1:n.diffexp, ])
  sd_no_tree <- matrixStats::colSds(no_tree_sim@count.matrix[1:n.diffexp, ])
  
  mean_us_tree <- colMeans(us_tree_sim@count.matrix[1:n.diffexp, ])
  sd_us_tree <-  matrixStats::colSds(us_tree_sim@count.matrix[1:n.diffexp, ])
  
  testthat::expect_equal(mean(mean_no_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_no_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(mean_us_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_us_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  testthat::expect_equal(mean(sd_us_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 2]),
                         mean(sd_us_tree[id_cond == 2]), tol = 1e-2)
  
})

test_that("Moments of the Negative Binomial and Log Normal are the same, with shift, random tree, with replicates", {
  
  tdir <- tempdir()
  
  ## Tree and conds
  n <- 20
  set.seed(1789)
  tree_norep <- ape::rphylo(n, 1, 0)
  
  r <- 2
  n <- r * n
  
  tree <- add_replicates(tree_norep, r)
  
  id_species <- sub("_.*$", "", tree$tip.label)
  id_species <- factor(id_species)
  names(id_species) <- tree$tip.label
  
  ## Condn alternate tips
  species_names <- unique(sub("_.*$", "", tree$tip.label))
  cond_species <- rep(c(1, 2), length(species_names) / 2)
  names(cond_species) <- species_names
  
  id_cond <- id_species
  id_cond <- cond_species[as.vector(id_cond)]
  id_cond <- as.factor(id_cond)
  names(id_cond) <- names(id_species)
  
  ## Parameters
  n.vars <- 5000
  relmean <- 100
  reldisp <- 0.1
  seqdepth <- 50000
  
  n.diffexp <- n.vars / 2
  effect.size <- 2
  
  true_mean_1 <- seqdepth / n.vars
  true_mean_2 <- 2 * seqdepth * (effect.size + 1) / n.vars / (effect.size + 2)
  true_var_1 <- reldisp * true_mean_1^2 + true_mean_1
  true_sd_1 <- sqrt(true_var_1)
  true_var_2 <- reldisp * true_mean_2^2 + true_mean_2
  true_sd_2 <- sqrt(true_var_2)
  
  ## No tree
  set.seed(18570823)
  dataset <- "no_tree"
  no_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE
  )
  
  ## us tree
  set.seed(18570823)
  dataset <- "tree"
  us_tree_sim <- phylocompcodeR::generateSyntheticData(dataset = dataset,
                                                       n.vars = n.vars,
                                                       samples.per.cond = n / 2,
                                                       n.diffexp = n.diffexp,
                                                       repl.id = 1,
                                                       seqdepth = seqdepth,
                                                       minfact = 1,
                                                       maxfact = 1,
                                                       relmeans = rep(relmean, n.vars),
                                                       dispersions = rep(reldisp, n.vars),
                                                       fraction.upregulated = 1,
                                                       between.group.diffdisp = FALSE,
                                                       filter.threshold.total = 1,
                                                       filter.threshold.mediancpm = 0,
                                                       fraction.non.overdispersed = 0,
                                                       random.outlier.high.prob = 0,
                                                       random.outlier.low.prob = 0,
                                                       single.outlier.high.prob = 0,
                                                       single.outlier.low.prob = 0,
                                                       effect.size = effect.size,
                                                       output.file = file.path(tdir, "tmp.rds"),
                                                       tree = tree,
                                                       id.condition = id_cond,
                                                       id.species = id_species,
                                                       lengths.relmeans = NULL,
                                                       lengths.dispersions = NULL,
                                                       lengths.phylo = FALSE,
                                                       prop.var.tree = 0.9
  )
  
  
  ## Comparisons
  mean_no_tree <- colMeans(no_tree_sim@count.matrix[1:n.diffexp, ])
  sd_no_tree <- matrixStats::colSds(no_tree_sim@count.matrix[1:n.diffexp, ])
  
  mean_us_tree <- colMeans(us_tree_sim@count.matrix[1:n.diffexp, ])
  sd_us_tree <-  matrixStats::colSds(us_tree_sim@count.matrix[1:n.diffexp, ])
  
  testthat::expect_equal(mean(mean_no_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_no_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(mean_us_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_us_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  testthat::expect_equal(mean(sd_us_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 2]),
                         mean(sd_us_tree[id_cond == 2]), tol = 1e-2)
  
})
