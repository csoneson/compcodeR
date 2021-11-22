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
  no_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  us_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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

  ## Correlation Comparisons
  corup <- function(x, ...) {
    cc <- cor(x, ...)
    return(cc[upper.tri(cc)])
  }
  # Pearson correlations are equal
  cor_no_tree <- corup(no_tree_sim@count.matrix)
  cor_tree <- corup(us_tree_sim@count.matrix)
  testthat::expect_equal(median(cor_no_tree), median(cor_tree), tol = 1e-2)
  # Spearman correlations are equal
  cor_no_tree <- corup(no_tree_sim@count.matrix, method = "spearman")
  cor_tree <- corup(us_tree_sim@count.matrix, method = "spearman")
  testthat::expect_equal(median(cor_no_tree), median(cor_tree), tol = 1e-2)

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
  no_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  us_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  
  true_no_tree <- no_tree_sim@variable.annotations$truemeans.S2[1:n.diffexp] / sum(no_tree_sim@variable.annotations$truemeans.S2) * seqdepth
  mean_real_no_tree <- mean(true_no_tree)
  sd_real_no_tree <- sd(true_no_tree)
  
  mean_us_tree <- colMeans(us_tree_sim@count.matrix[1:n.diffexp, ])
  sd_us_tree <-  matrixStats::colSds(us_tree_sim@count.matrix[1:n.diffexp, ])
  
  true_us_tree <- us_tree_sim@variable.annotations$truemeans.S2[1:n.diffexp] / sum(us_tree_sim@variable.annotations$truemeans.S2) * seqdepth
  mean_real_us_tree <- mean(true_us_tree)
  sd_real_us_tree <- sd(true_us_tree)
  
  testthat::expect_equal(mean(mean_no_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_no_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  testthat::expect_equal(mean(mean_no_tree[id_cond == 2]), mean_real_no_tree, tol = 1e-2)
  
  testthat::expect_equal(mean(mean_us_tree[id_cond == 1]), true_mean_1, tol = 1e-2)
  testthat::expect_equal(mean(mean_us_tree[id_cond == 2]), true_mean_2, tol = 1e-2)
  testthat::expect_equal(mean(mean_us_tree[id_cond == 2]), mean_real_no_tree, tol = 1e-2)
  
  testthat::expect_equal(mean(sd_no_tree[id_cond == 1]), true_sd_1, tol = 1e-2)
  testthat::expect_equal(mean(sd_us_tree[id_cond == 1]), true_sd_1, tol = 1e-2)

  testthat::expect_equal(mean(sd_no_tree[id_cond == 2]),
                         mean(sd_us_tree[id_cond == 2]), tol = 1e-2)
  
  ## Correlation Comparisons
  corup <- function(x, ...) {
    cc <- cor(x, ...)
    return(cc[upper.tri(cc)])
  }
  # Pearson correlations are equal
  cor_no_tree <- corup(no_tree_sim@count.matrix)
  cor_tree <- corup(us_tree_sim@count.matrix)
  testthat::expect_equal(median(cor_no_tree), median(cor_tree), tol = 1e-2)
  # Spearman correlations are equal
  cor_no_tree <- corup(no_tree_sim@count.matrix, method = "spearman")
  cor_tree <- corup(us_tree_sim@count.matrix, method = "spearman")
  testthat::expect_equal(median(cor_no_tree), median(cor_tree), tol = 1e-2)
  
})

test_that("Moments of the Negative Binomial and Log Normal are the same, varying mean", {
  
  tdir <- tempdir()
  
  ## Tree and conds
  n <- 14
  tree <- ape::stree(n)
  tree$edge.length <- rep(1, nrow(tree$edge))
  id_species <- tree$tip.label
  id_species <- factor(id_species)
  
  tree$tip.label <- paste0(tree$tip.label, "_rep")
  
  names(id_species) <- tree$tip.label
  
  ## Condn alternate tips
  species_names <- tree$tip.label
  cond_species <-  c(2, 1, 2, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1) #rep(c(1, 2), length(species_names) / 2)
  names(cond_species) <- species_names
  
  id_cond <- cond_species
  id_cond <- as.factor(id_cond)
  
  effect.size <- 0
  
  ## Parameters
  set.seed(18570823)
  n.vars <- 3560
  length.mu.phi.estimates <- readRDS(system.file("extdata", "Stern2018.Length.Mu.Phi.Estimates.rds", package = "compcodeR"))
  relmean <- length.mu.phi.estimates$stern2018.length.mu
  reldisp <- length.mu.phi.estimates$stern2018.length.phi
  seqdepth <- 2000000
  
  n.diffexp <- 0
  
  nrep <- 1
  
  ## No tree
  set.seed(18570821)
  dataset <- "no_tree"
  count_no_tree <- array(NA, dim = c(n.vars, n, nrep))
  for (rep in 1:nrep) {
    no_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
                                                    n.vars = n.vars,
                                                    samples.per.cond = n / 2,
                                                    n.diffexp = n.diffexp,
                                                    repl.id = 1,
                                                    seqdepth = seqdepth,
                                                    minfact = 1,
                                                    maxfact = 1,
                                                    relmeans = relmean,
                                                    dispersions = reldisp,
                                                    fraction.upregulated = 0.5,
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
                                                    # tree = tree,
                                                    # model.process = "BM",
                                                    id.condition = id_cond,
                                                    id.species = id_species,
                                                    lengths.relmeans = NULL,
                                                    lengths.dispersions = NULL,
                                                    lengths.phylo = FALSE
    )
    count_no_tree[, , rep] <- no_tree_sim@count.matrix
  }
  
  ## us tree
  set.seed(18570821)
  dataset <- "tree"
  count_tree <- array(NA, dim = c(n.vars, n, nrep))
  for (rep in 1:nrep) {
    us_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
                                                    n.vars = n.vars,
                                                    samples.per.cond = n / 2,
                                                    n.diffexp = n.diffexp,
                                                    repl.id = 1,
                                                    seqdepth = seqdepth,
                                                    minfact = 1,
                                                    maxfact = 1,
                                                    relmeans = relmean,
                                                    dispersions = reldisp,
                                                    fraction.upregulated = 0.5,
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
                                                    model.process = "BM",
                                                    prop.var.tree = 1.0,
                                                    id.condition = id_cond,
                                                    id.species = id_species,
                                                    lengths.relmeans = NULL,
                                                    lengths.dispersions = NULL,
                                                    lengths.phylo = FALSE
    )
    count_tree[, , rep] <- us_tree_sim@count.matrix
  }
  
  ## Correlation Comparisons
  corup <- function(x, ...) {
    cc <- cor(x, ...)
    return(cc[upper.tri(cc)])
  }
  # Pearson correlations are equal
  cor_no_tree <- rowMeans(apply(count_no_tree, 3, corup))
  cor_tree <- rowMeans(apply(count_tree, 3, corup))
  testthat::expect_equal(median(cor_no_tree), median(cor_tree), tol = 1e-2)
  # But Spearman are not correlations not are equal
  cor_no_tree <- rowMeans(apply(count_no_tree, 3, corup, method = "spearman"))
  cor_tree <- rowMeans(apply(count_tree, 3, corup, method = "spearman"))
  testthat::expect_true(median(cor_no_tree) <= median(cor_tree))
  
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
  no_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  us_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  no_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
  us_tree_sim <- compcodeR::generateSyntheticData(dataset = dataset,
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
