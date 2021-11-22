test_that("generateSyntheticData fails with wrong inputs", {
  tdir <- tempdir()
  
  expect_error(generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    output.file = file.path(tdir, "tmp.txt")
  ), "output.file must be an .rds file.")
  
  expect_error(generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    effect.size = 1:3, 
    output.file = file.path(tdir, "tmp.rds")
  ), "The length of the effect.size vector must be the same as the number of simulated genes.")
  
  expect_error(generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    relmeans = 1:3,
    output.file = file.path(tdir, "tmp.rds")
  ), "The length of the relmeans vector must be the same as the number of simulated genes.")
  
  expect_error(generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    dispersions = 1:3,
    output.file = file.path(tdir, "tmp.rds")
  ), "The number of provided dispersions must be the same as the number of simulated genes.")
})

test_that("compData object checks work", {
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 0, 
    output.file = NULL
  )
  
  expect_equal(checkDataObject(tmp), "Data object looks ok.")
  
  l <- convertcompDataToList(tmp)
  expect_is(l, "list")
  cpd <- convertListTocompData(l)
  expect_is(cpd, "compData")
  expect_equal(checkDataObject(cpd), "Data object looks ok.")
  l$count.matrix <- NULL
  expect_message({cpd2 <- convertListTocompData(l)}, 
                 "Cannot convert list to compData object")
  expect_equal(cpd2, NULL)
  
  tmp2 <- tmp; expect_error({count.matrix(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({sample.annotations(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({filtering(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({analysis.date(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({package.version(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({info.parameters(tmp2) <- 1:3})
  
  expect_equal(check_compData(tmp), TRUE)
  expect_equal(check_compData(count.matrix(tmp)), "This is not an S4 object.")
  tmp2 <- tmp; count.matrix(tmp2) <- as.matrix(numeric(0)); expect_equal(check_compData(tmp2), "Object must contain a non-empty count matrix.")
  tmp2 <- tmp; sample.annotations(tmp2) <- as.data.frame(NULL); expect_equal(check_compData(tmp2), "Object must contain a non-empty sample annotation data frame.")
  tmp2 <- tmp; sample.annotations(tmp2) <- data.frame(condition = numeric(0)); expect_equal(check_compData(tmp2), "The sample.annotations must contain a column named condition.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(); expect_equal(check_compData(tmp2), "Object must contain a non-empty list called info.parameters.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(tmp = 1); expect_equal(check_compData(tmp2), "info.parameters list must contain an entry named 'dataset'.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(dataset = ""); expect_equal(check_compData(tmp2), "info.parameters list must contain an entry named 'uID'.")
  tmp2 <- tmp; count.matrix(tmp2) <- count.matrix(tmp2)[1:10, ]; expect_equal(check_compData(tmp2), "count.matrix and variable.annotations do not contain the same number of rows.")
  tmp2 <- tmp; rownames(count.matrix(tmp2)) <- paste0("r", rownames(count.matrix(tmp2))); expect_equal(check_compData(tmp2), "The rownames of count.matrix and variable.annotations are not the same.")
  tmp2 <- tmp; count.matrix(tmp2) <- count.matrix(tmp2)[, 1:2]; expect_equal(check_compData(tmp2), "The number of columns of count.matrix is different from the number of rows of sample.annotations.")
  tmp2 <- tmp; colnames(count.matrix(tmp2)) <- paste0("r", colnames(count.matrix(tmp2))); expect_equal(check_compData(tmp2), "The colnames of count.matrix are different from the rownames of sample.annotations.")
  
  expect_equal(check_compData_results(tmp), "Object must contain a list named 'method.names' identifying the differential expression method used.")
})

test_that("phyloCompData object checks work", {
  tree <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  idsp <- as.factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idsp) <- tree$tip.label
  idcond <- c(1, 1, 1, 1, 2, 2, 2, 2)
  names(idcond) <- tree$tip.label
  
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 10,
    tree = tree,
    id.species =  idsp,
    id.condition = idcond,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = NULL
  )
  
  expect_equal(checkDataObject(tmp), "Data object looks ok.")
  
  l <- convertphyloCompDataToList(tmp)
  expect_is(l, "list")

  cpd <- convertListTocompData(l)
  expect_is(cpd, "compData")
  expect_equal(checkDataObject(cpd), "Data object looks ok.")
  expect_null(length.matrix(cpd))
  expect_null(phylo.tree(cpd))
  expect_error(phylo.tree(cpd) <- l$tree, "There is no 'phylo.tree' slot in a 'compData' object. Please use a 'phyloCompData' object.")
  expect_error(length.matrix(cpd) <- l$length.matrix, "There is no 'lenght.matrix' slot in a 'compData' object. Please use a 'phyloCompData' object.")
  
  
  cpd <- convertListTophyloCompData(l)
  expect_is(cpd, "phyloCompData")
  expect_equal(checkDataObject(cpd), "Data object looks ok.")
  
  cpdbis <- phyloCompData(l$count.matrix, l$sample.annotations, 
                          l$info.parameters, l$variable.annotations, 
                          l$filtering, character(0), 
                          l$package.version, l$method.names, 
                          l$code, l$result.table,
                          l$tree,
                          l$length.matrix)
  expect_equal(cpd, cpdbis)
  
  l$count.matrix <- NULL
  expect_message({cpd2 <- convertListTocompData(l)}, 
                 "Cannot convert list to compData object")
  expect_message({cpd2 <- convertListTophyloCompData(l)}, 
                 "Cannot convert list to compData object")
  expect_equal(cpd2, NULL)
  
  tmp2 <- tmp; expect_error({count.matrix(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({sample.annotations(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({filtering(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({analysis.date(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({package.version(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({info.parameters(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({length.matrix(tmp2) <- 1:3})
  tmp2 <- tmp; expect_error({phylo.tree(tmp2) <- 1:3})
  
  expect_equal(check_compData(tmp), TRUE)
  expect_equal(check_compData(count.matrix(tmp)), "This is not an S4 object.")
  tmp2 <- tmp; count.matrix(tmp2) <- as.matrix(numeric(0)); expect_equal(check_phyloCompData(tmp2), "Object must contain a non-empty count matrix.")
  tmp2 <- tmp; sample.annotations(tmp2) <- as.data.frame(NULL); expect_equal(check_phyloCompData(tmp2), "Object must contain a non-empty sample annotation data frame.")
  tmp2 <- tmp; sample.annotations(tmp2) <- data.frame(condition = numeric(0)); expect_equal(check_phyloCompData(tmp2), "The sample.annotations must contain a column named condition.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(); expect_equal(check_phyloCompData(tmp2), "Object must contain a non-empty list called info.parameters.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(tmp = 1); expect_equal(check_phyloCompData(tmp2), "info.parameters list must contain an entry named 'dataset'.")
  tmp2 <- tmp; info.parameters(tmp2) <- list(dataset = ""); expect_equal(check_phyloCompData(tmp2), "info.parameters list must contain an entry named 'uID'.")
  tmp2 <- tmp; count.matrix(tmp2) <- count.matrix(tmp2)[1:10, ]; expect_equal(check_phyloCompData(tmp2), "count.matrix and variable.annotations do not contain the same number of rows.")
  tmp2 <- tmp; rownames(count.matrix(tmp2)) <- paste0("r", rownames(count.matrix(tmp2))); expect_equal(check_phyloCompData(tmp2), "The rownames of count.matrix and variable.annotations are not the same.")
  tmp2 <- tmp; count.matrix(tmp2) <- count.matrix(tmp2)[, 1:2]; expect_equal(check_phyloCompData(tmp2), "The number of columns of count.matrix is different from the number of rows of sample.annotations.")
  tmp2 <- tmp; colnames(count.matrix(tmp2)) <- paste0("r", colnames(count.matrix(tmp2))); expect_equal(check_phyloCompData(tmp2), "The colnames of count.matrix are different from the rownames of sample.annotations.")
  
  tmp2 <- tmp; length.matrix(tmp2) <- length.matrix(tmp2)[1:10, ]; expect_equal(check_phyloCompData(tmp2), "The dimension of count.matrix is different from the dimension of length.matrix.")
  tmp2 <- tmp; rownames(length.matrix(tmp2)) <- paste0("r", rownames(length.matrix(tmp2))); expect_equal(check_phyloCompData(tmp2), "The rownames of count.matrix are different from the rownames of length.matrix.")
  tmp2 <- tmp; length.matrix(tmp2) <- length.matrix(tmp2)[, 1:2]; expect_equal(check_phyloCompData(tmp2), "The dimension of count.matrix is different from the dimension of length.matrix.")
  tmp2 <- tmp; colnames(length.matrix(tmp2)) <- paste0("r", colnames(length.matrix(tmp2))); expect_equal(check_phyloCompData(tmp2), "The colnames of count.matrix are different from the colnames of length.matrix.")
  
  tmp2 <- tmp; phylo.tree(tmp2)$tip.label <- NULL; expect_equal(check_phyloCompData(tmp2), "The tips of the phylogeny are not named.")
  tmp2 <- tmp; phylo.tree(tmp2)$tip.label <- paste0("r", phylo.tree(tmp2)$tip.label); expect_equal(check_phyloCompData(tmp2), "Column names of count.matrix do not match the tip labels.")
  tmp2 <- tmp; phylo.tree(tmp2)$tip.label <- phylo.tree(tmp2)$tip.label[c(2, 1, 3:8)]; expect_equal(check_phyloCompData(tmp2), "Column names of count.matrix do not match the tip labels.")
  tmp2 <- tmp; colnames(count.matrix(tmp2)) <- colnames(count.matrix(tmp2))[c(2, 1, 3:8)]; expect_equal(check_phyloCompData(tmp2), "The colnames of count.matrix are different from the rownames of sample.annotations.")
  tmp2 <- tmp; colnames(length.matrix(tmp2)) <- colnames(length.matrix(tmp2))[c(2, 1, 3:8)]; expect_equal(check_phyloCompData(tmp2), "The colnames of count.matrix are different from the colnames of length.matrix.")
  tmp2 <- tmp; rownames(sample.annotations(tmp2)) <- rownames(sample.annotations(tmp2))[c(2, 1, 3:8)]; expect_equal(check_phyloCompData(tmp2), "The colnames of count.matrix are different from the rownames of sample.annotations.")
  tmp2 <- tmp; sample.annotations(tmp2)$id.species <- rep(1, 8); expect_equal(check_phyloCompData(tmp2), "Error in checkSpecies(ids, \"id.species\", phylo.tree(object), tol = 1e-10,  : \n  The provided species do not match with the tree branch lengths. Please check the 'id.species' vector.\n")
  tmp2 <- tmp; sample.annotations(tmp2)$id.species <- NULL; expect_equal(check_phyloCompData(tmp2), "The sample.annotations must contain a column named id.species.")
  
  expect_equal(check_compData_results(tmp), "Object must contain a list named 'method.names' identifying the differential expression method used.")
})

test_that("generateSyntheticData works", {
  tdir <- tempdir()
  
  ## No DEGs
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 0, 
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(variable.annotations(tmp)$truedispersions.S1, 
               variable.annotations(tmp)$truedispersions.S2)
  expect_equal(variable.annotations(tmp)$truemeans.S1,
               variable.annotations(tmp)$truemeans.S2)
  expect_equal(sum(variable.annotations(tmp)$n.random.outliers.up.S1 + 
                     variable.annotations(tmp)$n.random.outliers.up.S2 + 
                     variable.annotations(tmp)$n.random.outliers.down.S1 + 
                     variable.annotations(tmp)$n.random.outliers.down.S2 + 
                     variable.annotations(tmp)$n.single.outliers.up.S1 + 
                     variable.annotations(tmp)$n.single.outliers.up.S2 + 
                     variable.annotations(tmp)$n.single.outliers.down.S1 + 
                     variable.annotations(tmp)$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(variable.annotations(tmp)$truelog2foldchanges)), 0)
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation + 
                     variable.annotations(tmp)$differential.expression), 0)
               
  ## Specify effect sizes individually
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    effect.size = c(exp(abs(rnorm(5))), exp(-abs(rnorm(5))), rep(1, 40)),
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(variable.annotations(tmp)$truedispersions.S1, 
               variable.annotations(tmp)$truedispersions.S2)
  expect_equal(variable.annotations(tmp)$truemeans.S1[-(1:10)],
               variable.annotations(tmp)$truemeans.S2[-(1:10)])
  expect_equal(sum(variable.annotations(tmp)$n.random.outliers.up.S1 + 
                     variable.annotations(tmp)$n.random.outliers.up.S2 + 
                     variable.annotations(tmp)$n.random.outliers.down.S1 + 
                     variable.annotations(tmp)$n.random.outliers.down.S2 + 
                     variable.annotations(tmp)$n.single.outliers.up.S1 + 
                     variable.annotations(tmp)$n.single.outliers.up.S2 + 
                     variable.annotations(tmp)$n.single.outliers.down.S1 + 
                     variable.annotations(tmp)$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(variable.annotations(tmp)$truelog2foldchanges[-(1:10)])), 0)
  expect_equal(sign(variable.annotations(tmp)$truelog2foldchanges),
               c(rep(1, 5), rep(-1, 5), rep(0, 40)))
  expect_equal(variable.annotations(tmp)$upregulation, 
               c(rep(1, 5), rep(0, 45)))
  expect_equal(variable.annotations(tmp)$downregulation, 
               c(rep(0, 5), rep(1, 5), rep(0, 40)))
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation), 10)
  
  ## Different dispersions between groups
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    between.group.diffdisp = TRUE,
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(any(variable.annotations(tmp)$truedispersions.S1 !=  
                     variable.annotations(tmp)$truedispersions.S2),
               TRUE)
  expect_equal(variable.annotations(tmp)$upregulation, 
               c(rep(1, 10), rep(0, 40)))
  expect_equal(variable.annotations(tmp)$downregulation, 
               rep(0, 50))
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation), 10)
  
  ## Not overdispersed
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    between.group.diffdisp = FALSE, 
    fraction.non.overdispersed = 0.5,
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(variable.annotations(tmp)$truedispersions.S1,  
               variable.annotations(tmp)$truedispersions.S2)
  expect_equal(any(variable.annotations(tmp)$truedispersions.S1 == 0), TRUE)
  expect_equal(info.parameters(tmp)$fraction.non.overdispersed, 0.5)

  ## Outliers
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    random.outlier.high.prob = 0.1,
    random.outlier.low.prob = 0.1,
    single.outlier.high.prob = 0.1,
    single.outlier.low.prob = 0.1,
    output.file = NULL
  )
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.up.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.up.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.down.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.down.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.up.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.up.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.down.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.down.S2 > 0), TRUE)
  
  ## Summary report
  expect_error(summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp.rds")),
               "output.file must be an .html file.")
  summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp_summaryrep.html"))
  expect_equal(file.exists(normalizePath(file.path(tdir, "tmp_summaryrep.html"), 
                                         winslash = "/")), TRUE)
})

test_that("generateSyntheticData works - with lengths and phylo", {
  tdir <- tempdir()
  
  ## Errors with lengths
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 50, 
      samples.per.cond = 4, n.diffexp = 5, 
      repl.id = 1, 
      lengths.relmeans = rpois(40, 1e4),
      lengths.dispersions = rgamma(50, 1, 1),
      output.file = NULL
    ),
    "The length of the 'lengths.relmeans' vector must be the same as the number of simulated genes.")
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 50, 
      samples.per.cond = 4, n.diffexp = 5, 
      repl.id = 1, 
      lengths.relmeans = rpois(50, 1e4),
      lengths.dispersions = rgamma(40, 1, 1),
      output.file = NULL
    ),
    "The length of the 'lengths.dispersions' vector must be the same as the number of simulated genes.")
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      lengths.relmeans = rpois(50, 1e4),
      output.file = NULL
    ),
    "For lengths to be used, both the 'lengths.relmeans' and 'lengths.dispersions' vectors must be provided.")
  
  ## Errors and warnings with tree
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1, 
      tree = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);",
      output.file = NULL
    ),
    "The `tree` must be of class `phylo` from package `ape`.")
  
  tree <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):0.5,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      output.file = NULL
    ),
    "The tree should be ultrametric.")
  
  tree <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 5, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      output.file = NULL
    ),
    "The tree should have as many species as `samples.per.cond` times two")
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      id.species =  as.factor(c("A", "A", "A", "B", "C", "C", "D")),
      output.file = NULL
    ),
    "`id.species` should have the same length as the number of taxa in the tree.")
  expect_warning(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      id.species =  c("A", "A", "A", "B", "C", "C", "D", "D"),
      output.file = NULL
    ),
    "Vector 'id.species' must be a factor.")
  expect_warning(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      id.species =  as.factor(c("A", "A", "A", "B", "C", "C", "D", "D")),
      output.file = NULL
    ),
    "`id.species` is not named. I'm naming them, assuming they are in the same order as the tree.")
  
  idsp <- as.factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idsp) <- c("F", tree$tip.label[-1])
  
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      id.species =  idsp,
      output.file = NULL
    ),
    "`id.species` names do not match the tip labels.")
  
  idsp <- as.factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idsp) <- tree$tip.label
  
  expect_warning(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 500, 
      samples.per.cond = 4, n.diffexp = 50, 
      repl.id = 1,
      tree = tree,
      id.condition = c(1, 1, 1, 1, 2, 2, 2, 2),
      id.species =  idsp,
      output.file = NULL
    ),
    "`id.condition` is not named. I'm naming them, assuming they are in the same order as the tree.")
  
  idcond <- c(1, 1, 1, 1, 2, 2, 2, 2)
  names(idcond) <- tree$tip.label
  
  ## No DEGs
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 0,
    tree = tree,
    id.species =  idsp,
    id.condition = idcond,
    output.file = NULL
  )
  expect_is(tmp, "phyloCompData")
  expect_equal(variable.annotations(tmp)$truedispersions.S1, 
               variable.annotations(tmp)$truedispersions.S2)
  expect_equal(variable.annotations(tmp)$truemeans.S1,
               variable.annotations(tmp)$truemeans.S2)
  expect_equal(sum(variable.annotations(tmp)$n.random.outliers.up.S1 + 
                     variable.annotations(tmp)$n.random.outliers.up.S2 + 
                     variable.annotations(tmp)$n.random.outliers.down.S1 + 
                     variable.annotations(tmp)$n.random.outliers.down.S2 + 
                     variable.annotations(tmp)$n.single.outliers.up.S1 + 
                     variable.annotations(tmp)$n.single.outliers.up.S2 + 
                     variable.annotations(tmp)$n.single.outliers.down.S1 + 
                     variable.annotations(tmp)$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(variable.annotations(tmp)$truelog2foldchanges)), 0)
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation + 
                     variable.annotations(tmp)$differential.expression), 0)
  
  ## Specify effect sizes individually
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 10,
    tree = tree,
    id.species =  idsp,
    id.condition = idcond,
    effect.size = c(exp(abs(rnorm(5))), exp(-abs(rnorm(5))), rep(1, 40)),
    output.file = NULL
  )
  expect_is(tmp, "phyloCompData")
  expect_equal(variable.annotations(tmp)$truedispersions.S1, 
               variable.annotations(tmp)$truedispersions.S2)
  expect_equal(variable.annotations(tmp)$truemeans.S1[-(1:10)],
               variable.annotations(tmp)$truemeans.S2[-(1:10)])
  expect_equal(sum(variable.annotations(tmp)$n.random.outliers.up.S1 + 
                     variable.annotations(tmp)$n.random.outliers.up.S2 + 
                     variable.annotations(tmp)$n.random.outliers.down.S1 + 
                     variable.annotations(tmp)$n.random.outliers.down.S2 + 
                     variable.annotations(tmp)$n.single.outliers.up.S1 + 
                     variable.annotations(tmp)$n.single.outliers.up.S2 + 
                     variable.annotations(tmp)$n.single.outliers.down.S1 + 
                     variable.annotations(tmp)$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(variable.annotations(tmp)$truelog2foldchanges[-(1:10)])), 0)
  expect_equal(sign(variable.annotations(tmp)$truelog2foldchanges),
               c(rep(1, 5), rep(-1, 5), rep(0, 40)))
  expect_equal(variable.annotations(tmp)$upregulation, 
               c(rep(1, 5), rep(0, 45)))
  expect_equal(variable.annotations(tmp)$downregulation, 
               c(rep(0, 5), rep(1, 5), rep(0, 40)))
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation), 10)
  
  ## Different dispersions between groups
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 10,
    between.group.diffdisp = TRUE,
    tree = tree,
    id.species =  idsp,
    id.condition = idcond,
    output.file = NULL
  )
  expect_is(tmp, "phyloCompData")
  expect_equal(any(variable.annotations(tmp)$truedispersions.S1 !=  
                     variable.annotations(tmp)$truedispersions.S2),
               TRUE)
  expect_equal(variable.annotations(tmp)$upregulation, 
               c(rep(1, 10), rep(0, 40)))
  expect_equal(variable.annotations(tmp)$downregulation, 
               rep(0, 50))
  expect_equal(sum(variable.annotations(tmp)$upregulation + 
                     variable.annotations(tmp)$downregulation), 10)
  
  ## Not overdispersed
  expect_error(
    generateSyntheticData(
      dataset = "B_625_625", n.vars = 50, 
      samples.per.cond = 4, n.diffexp = 10,
      between.group.diffdisp = FALSE, 
      fraction.non.overdispersed = 0.5,
      tree = tree,
      id.species =  idsp,
      id.condition = idcond,
      output.file = NULL
    ),
    "The Phylogenetic Poisson lognormal distribution is always over-dispersed.")
  
  ## Outliers
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 10,
    tree = tree,
    id.species =  idsp,
    random.outlier.high.prob = 0.1,
    random.outlier.low.prob = 0.1,
    single.outlier.high.prob = 0.1,
    single.outlier.low.prob = 0.1,
    output.file = NULL
  )
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.up.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.up.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.down.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.random.outliers.down.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.up.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.up.S2 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.down.S1 > 0), TRUE)
  expect_equal(any(variable.annotations(tmp)$n.single.outliers.down.S2 > 0), TRUE)
  
  ## Summary report
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 4, n.diffexp = 10,
    tree = tree,
    id.species =  idsp,
    id.condition = idcond,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = NULL
  )
  
  expect_error(summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp.rds")),
               "output.file must be an .html file.")
  summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp_summaryrep.html"))
  expect_equal(file.exists(normalizePath(file.path(tdir, "tmp_summaryrep.html"), 
                                         winslash = "/")), TRUE)
})

test_that("help functions work", {
  listcreateRmd()
  
  expect_warning(expect_error(checkRange("hello", "name", 0, 1), "Illegal value"), "NAs introduced by coercion")
  expect_equal(checkRange(-1, "name", 0, 1), 0)
  expect_equal(checkRange(2, "name", 0, 1), 1)
  expect_equal(checkRange("-1", "name", 0, 1), 0)
  
  expect_equal(shorten.method.names(c("AUC", "ROC, all replicates")),
               c("auc", "rocall"))
  
  set.seed(1)
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    output.file = NULL
  )
  
  expect_error(
    computeMval(count.matrix(tmp), c(3, sample.annotations(tmp)$condition[-1])),
    "Must have exactly two groups to calculate M-value"
  )
  expect_error(
    computeAval(count.matrix(tmp), c(3, sample.annotations(tmp)$condition[-1])),
    "Must have exactly two groups to calculate A-value"
  )
  
  mval <- computeMval(count.matrix(tmp), sample.annotations(tmp)$condition)
  aval <- computeAval(count.matrix(tmp), sample.annotations(tmp)$condition)
  
  expect_is(mval, "numeric")
  expect_is(aval, "numeric")
  expect_equal(length(mval), nrow(count.matrix(tmp)))
  expect_equal(length(aval), nrow(count.matrix(tmp)))
})

test_that("runDiffExp works", {
  
  tdir <- tempdir()
  set.seed(1)  ## note that with other seeds, the number of genes 
               ## passing the filtering threshold could be different
  
  testdat <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    output.file = file.path(tdir, "B_625_625_5spc_repl1.rds")
  )
  
  expect_equal(checkDataObject(testdat), "Data object looks ok.")
  
  tmp <- readRDS(normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"))
  expect_is(tmp, "compData")
  expect_is(count.matrix(tmp), "matrix")
  expect_equal(nrow(count.matrix(tmp)), 500)
  expect_equal(ncol(count.matrix(tmp)), 10)
  expect_is(sample.annotations(tmp), "data.frame")
  expect_equal(nrow(sample.annotations(tmp)), 10)
  expect_equal(ncol(sample.annotations(tmp)), 2)
  expect_is(variable.annotations(tmp), "data.frame")
  expect_equal(nrow(variable.annotations(tmp)), 500)
  expect_equal(ncol(variable.annotations(tmp)), 18)
  expect_is(info.parameters(tmp), "list")
  expect_equal(info.parameters(tmp)$n.diffexp, 50)
  expect_equal(info.parameters(tmp)$dataset, "B_625_625")
  expect_equal(info.parameters(tmp)$fraction.upregulated, 0.5)
  expect_equal(info.parameters(tmp)$between.group.diffdisp, FALSE)
  expect_equal(info.parameters(tmp)$filter.threshold.total, 1)
  expect_equal(info.parameters(tmp)$filter.threshold.mediancpm, 0)
  expect_equal(info.parameters(tmp)$fraction.non.overdispersed, 0)
  expect_equal(info.parameters(tmp)$random.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$random.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$effect.size, 1.5)
  expect_equal(info.parameters(tmp)$samples.per.cond, 5)
  expect_equal(info.parameters(tmp)$repl.id, 1)
  expect_equal(info.parameters(tmp)$seqdepth, 1e5)
  expect_equal(info.parameters(tmp)$minfact, 0.7)
  expect_equal(info.parameters(tmp)$maxfact, 1.4)
  expect_equal(filtering(tmp), "total count >= 1 ;  median cpm >= 0")
  expect_is(analysis.date(tmp), "character")
  expect_equal(analysis.date(tmp), "")
  expect_is(package.version(tmp), "character")
  expect_equal(package.version(tmp), "")
  expect_is(method.names(tmp), "list")
  expect_equal(method.names(tmp), list())
  
  if (requireNamespace("baySeq", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "baySeq",
      Rmdfunction = "baySeq.createRmd",
      output.directory = tdir, norm.method = "edgeR",
      equaldisp = TRUE
    )
  }
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "DESeq2",
      Rmdfunction = "DESeq2.createRmd",
      output.directory = tdir, fit.type = "parametric",
      test = "Wald", beta.prior = TRUE,
      independent.filtering = TRUE, cooks.cutoff = TRUE,
      impute.outliers = TRUE
    )
  }
  if (requireNamespace("DSS", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "DSS",
      Rmdfunction = "DSS.createRmd",
      output.directory = tdir, norm.method = "quantile",
      disp.trend = FALSE
    )
  }
  if (requireNamespace("EBSeq", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "EBSeq",
      Rmdfunction = "EBSeq.createRmd",
      output.directory = tdir, norm.method = "median"
    )
  }
  # edgeR is in Imports -> always installed
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
    result.extent = "edgeR.exact",
    Rmdfunction = "edgeR.exact.createRmd",
    output.directory = tdir, norm.method = "TMM",
    trend.method = "movingave", disp.type = "tagwise"
  )
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
    result.extent = "edgeR.GLM",
    Rmdfunction = "edgeR.GLM.createRmd",
    output.directory = tdir, norm.method = "TMM",
    disp.type = "tagwise", disp.method = "CoxReid",
    trended = TRUE
  )
  # limma is in Imports -> always installed
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
    result.extent = "logcpm.limma",
    Rmdfunction = "logcpm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  if (requireNamespace("NBPSeq", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "NBPSeq",
      Rmdfunction = "NBPSeq.createRmd",
      output.directory = tdir, norm.method = "TMM",
      disp.method = "NBP"
    )
  }
  if (requireNamespace("NOISeq", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "NOISeq",
      Rmdfunction = "NOISeq.prenorm.createRmd",
      output.directory = tdir, norm.method = "TMM"
    )
  }
  # limma is in Imports -> always installed
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
    result.extent = "sqrtcpm.limma",
    Rmdfunction = "sqrtcpm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  if (requireNamespace("TCC", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "TCC",
      Rmdfunction = "TCC.createRmd",
      output.directory = tdir, norm.method = "tmm",
      test.method = "edger", iteration = 3,
      normFDR = 0.1, floorPDEG = 0.05
    )
  }
  if (requireNamespace("genefilter", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "ttest",
      Rmdfunction = "ttest.createRmd",
      output.directory = tdir, norm.method = "TMM"
    )
  }
  # limma is in Imports -> always installed
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "voom.limma",
    Rmdfunction = "voom.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  if (requireNamespace("genefilter", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
      result.extent = "voom.ttest",
      Rmdfunction = "voom.ttest.createRmd",
      output.directory = tdir, norm.method = "TMM"
    )
  }
  
  methods <- c("baySeq", "DESeq2", "DSS", "EBSeq", "edgeR.exact",
               "edgeR.GLM", "logcpm.limma", "NBPSeq", "NOISeq", 
               "sqrtcpm.limma", "TCC", "ttest", "voom.limma",
               "voom.ttest")
  pkgs <- c("baySeq", "DESeq2", "DSS", "EBSeq", "edgeR",
            "edgeR", "limma", "NBPSeq", "NOISeq", 
            "limma", "TCC", "genefilter", "limma",
            "genefilter")
  
  ## Test show() method
  m <- "edgeR.exact" # edgeR always installed
  tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")),
                               winslash = "/"))
  show(tmp)
  count.matrix(tmp) <- count.matrix(tmp)[, 1:4]
  show(tmp)
  
  for (i in seq_len(length(methods))) {
    m <- methods[i]
    pkg <- pkgs[i]
    if (requireNamespace(pkg, quietly = TRUE)) {
      tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), winslash = "/"))
      
      expect_is(tmp, "compData")
      expect_is(result.table(tmp), "data.frame")
      expect_equal(nrow(result.table(tmp)), 500)
      expect_is(code(tmp), "character")
      expect_is(analysis.date(tmp), "character")
      expect_is(compcodeR:::package.version(tmp), "character")
      expect_is(method.names(tmp), "list")
      expect_named(method.names(tmp), c("short.name", "full.name"))
      
      tmp2 <- tmp; result.table(tmp2) <- result.table(tmp2)[1:10, ]; expect_equal(check_compData_results(tmp2), "result.table must have the same number of rows as count.matrix.")
      tmp2 <- tmp; result.table(tmp2) <- data.frame(); expect_equal(check_compData_results(tmp2), "Object must contain a data frame named 'result.table'.")
      tmp2 <- tmp; result.table(tmp2)$score <- NULL; expect_equal(check_compData_results(tmp2), "result.table must contain a column named 'score'.")
    }
  }
  
  for (i in seq_len(length(methods))) {
    m <- methods[i]
    pkg <- pkgs[i]
    if (requireNamespace(pkg, quietly = TRUE)) {
      generateCodeHTMLs(
        normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), 
                      winslash = "/"), normalizePath(tdir)
      )
      expect_true(file.exists(normalizePath(file.path(
        tdir, paste0("B_625_625_5spc_repl1_", 
                     m, "_code.html")), winslash = "/")))
    }
  }
  
  ## Comparison report
  file.table <- data.frame(input.files = normalizePath(file.path(
    tdir, paste0("B_625_625_5spc_repl1_", 
                 c("voom.limma", "sqrtcpm.limma", "edgeR.exact", "edgeR.GLM"), ## Only packages in import
                 ".rds")), winslash = "/"))
  parameters <- NULL
  comp <- runComparison(file.table = file.table, output.directory = tdir,
                        parameters = parameters)
  
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = parameters, save.result.table = FALSE, knit.results = FALSE),
               "At least on of 'save.result.table' or 'knit.results' must be set to TRUE, otherwise the function does not produce anything.")
  comp2 <- runComparison(file.table = file.table, output.directory = tdir,
                         parameters = parameters, save.result.table = TRUE, knit.results = FALSE)
  ff <- list.files(path = tdir, pattern = "compcodeR_result_table_.*.rds", full.names = TRUE)
  resTable <- readRDS(ff[length(ff)])
  expect_equal(resTable$fp + resTable$tp + resTable$fn + resTable$tn, rep(500, nrow(resTable)))
  expect_equal(ncol(resTable), 14)
  expect_equal(nrow(resTable), 4)
  
  parameters <- list()
  par2 <- parameters; par2$incl.dataset <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with datasets")
  par2 <- parameters; par2$incl.nbr.samples <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with nbr.samples")
  par2 <- parameters; par2$incl.replicates <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with replicates")
  par2 <- parameters; par2$incl.de.methods <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with DE methods")
  
})

test_that("runDiffExp works - with lengths", {
  
  tdir <- tempdir()
  set.seed(1)  ## note that with other seeds, the number of genes 
               ## passing the filtering threshold could be different
  
  testdat <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 5, n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    id.species = factor(1:10),
    lengths.relmeans = rpois(500, 1e4),
    lengths.dispersions = rgamma(500, 1, 1),
    output.file = file.path(tdir, "B_625_625_5spc_repl1.rds")
  )
  
  expect_equal(checkDataObject(testdat), "Data object looks ok.")
  
  tmp <- readRDS(normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"))
  expect_is(tmp, "compData")
  expect_is(count.matrix(tmp), "matrix")
  expect_equal(nrow(count.matrix(tmp)), 498) # some are filtered
  expect_equal(ncol(count.matrix(tmp)), 10)
  expect_is(sample.annotations(tmp), "data.frame")
  expect_equal(nrow(sample.annotations(tmp)), 10)
  expect_equal(ncol(sample.annotations(tmp)), 2)
  expect_is(variable.annotations(tmp), "data.frame")
  expect_equal(nrow(variable.annotations(tmp)), 498)
  expect_equal(ncol(variable.annotations(tmp)), 22)
  expect_is(info.parameters(tmp), "list")
  expect_equal(info.parameters(tmp)$n.diffexp, 50)
  expect_equal(info.parameters(tmp)$dataset, "B_625_625")
  expect_equal(info.parameters(tmp)$fraction.upregulated, 0.5)
  expect_equal(info.parameters(tmp)$between.group.diffdisp, FALSE)
  expect_equal(info.parameters(tmp)$filter.threshold.total, 1)
  expect_equal(info.parameters(tmp)$filter.threshold.mediancpm, 0)
  expect_equal(info.parameters(tmp)$fraction.non.overdispersed, 0)
  expect_equal(info.parameters(tmp)$random.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$random.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$effect.size, 1.5)
  expect_equal(info.parameters(tmp)$samples.per.cond, 5)
  expect_equal(info.parameters(tmp)$repl.id, 1)
  expect_equal(info.parameters(tmp)$seqdepth, 1e5)
  expect_equal(info.parameters(tmp)$minfact, 0.7)
  expect_equal(info.parameters(tmp)$maxfact, 1.4)
  expect_equal(info.parameters(tmp)$nEff, 2.5)
  expect_equal(info.parameters(tmp)$nEffRatio, 1.0)
  expect_equal(filtering(tmp), "total count >= 1 ;  median cpm >= 0")
  expect_is(analysis.date(tmp), "character")
  expect_equal(analysis.date(tmp), "")
  expect_is(package.version(tmp), "character")
  expect_equal(package.version(tmp), "")
  expect_is(method.names(tmp), "list")
  expect_equal(method.names(tmp), list())
  
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    expect_message(
      runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "DESeq2.length",
      Rmdfunction = "DESeq2.length.createRmd",
      output.directory = tdir, fit.type = "parametric",
      test = "Wald", beta.prior = TRUE,
      independent.filtering = TRUE, cooks.cutoff = TRUE,
      impute.outliers = TRUE,
      extra.design.covariates = NULL,
      nas_as_ones = FALSE
    ),
    "there might be some NAs in the adjusted p values computed by DESeq2"
    )
  }
  # limma is in Imports
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "lengthNorm.limma",
    Rmdfunction = "lengthNorm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM",
    extra.design.covariates = NULL,
    length.normalization = "RPKM",
    data.transformation = "log2",
    block.factor = NULL
  )
  # phylolm is in Imports
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "phylolm",
    Rmdfunction = "phylolm.createRmd",
    output.directory = tdir, norm.method = "TMM",
    model = "BM", measurement_error = TRUE,
    extra.design.covariates = NULL,
    length.normalization = "RPKM",
    data.transformation = "log2"
  )
  
  methods <- c("DESeq2.length", "lengthNorm.limma", "phylolm")
  
  ## Test show() method
  m <- "lengthNorm.limma"
  tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")),
                               winslash = "/"))
  show(tmp)
  count.matrix(tmp) <- count.matrix(tmp)[, 1:4]
  show(tmp)
  
  for (m in methods) {
    if (m != "DESeq2.length" || requireNamespace("DESeq2", quietly = TRUE)) {
      tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), winslash = "/"))
      
      expect_is(tmp, "phyloCompData")
      expect_is(result.table(tmp), "data.frame")
      expect_equal(nrow(result.table(tmp)), 498)
      expect_is(code(tmp), "character")
      expect_is(analysis.date(tmp), "character")
      expect_is(compcodeR:::package.version(tmp), "character")
      expect_is(method.names(tmp), "list")
      expect_named(method.names(tmp), c("short.name", "full.name"))
      
      tmp2 <- tmp; result.table(tmp2) <- result.table(tmp2)[1:10, ]; expect_equal(check_compData_results(tmp2), "result.table must have the same number of rows as count.matrix.")
      tmp2 <- tmp; result.table(tmp2) <- data.frame(); expect_equal(check_compData_results(tmp2), "Object must contain a data frame named 'result.table'.")
      tmp2 <- tmp; result.table(tmp2)$score <- NULL; expect_equal(check_compData_results(tmp2), "result.table must contain a column named 'score'.")
    }
  }
  
  for (m in methods) {
    if (m != "DESeq2.length" || requireNamespace("DESeq2", quietly = TRUE)) {
      generateCodeHTMLs(
        normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), 
                      winslash = "/"), normalizePath(tdir)
      )
      expect_true(file.exists(normalizePath(file.path(
        tdir, paste0("B_625_625_5spc_repl1_", 
                     m, "_code.html")), winslash = "/")))
    }
  }
  
  ## Comparison report
  file.table <- data.frame(input.files = normalizePath(file.path(
    tdir, paste0("B_625_625_5spc_repl1_", 
                 methods[-1],
                 ".rds")), winslash = "/"))
  parameters <- NULL
  comp <- runComparison(file.table = file.table, output.directory = tdir,
                        parameters = parameters)
  
  comp2 <- runComparison(file.table = file.table, output.directory = tdir,
                         parameters = parameters, save.result.table = TRUE, knit.results = FALSE)
  ff <- list.files(path = tdir, pattern = "compcodeR_result_table_.*.rds", full.names = TRUE)
  resTable <- readRDS(ff[length(ff)])
  expect_equal(resTable$fp + resTable$tp + resTable$fn + resTable$tn, rep(498, nrow(resTable)))
  expect_equal(ncol(resTable), 14)
  expect_equal(nrow(resTable), length(methods) - 1)
  
  parameters <- list()
  par2 <- parameters; par2$incl.dataset <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with datasets")
  par2 <- parameters; par2$incl.nbr.samples <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with nbr.samples")
  par2 <- parameters; par2$incl.replicates <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with replicates")
  par2 <- parameters; par2$incl.de.methods <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with DE methods")
})

test_that("runDiffExp works - phylo", {
  
  tdir <- tempdir()
  
  tree <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")

  idsp <- as.factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idsp) <- tree$tip.label
  
  idcond <- c(1, 1, 1, 1, 2, 2, 2, 2)
  names(idcond) <- tree$tip.label
  
  set.seed(1)  ## note that with other seeds, the number of genes 
               ## passing the filtering threshold could be different
  
  testdat <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 500, 
    samples.per.cond = 4, n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    tree = tree,
    id.condition = idcond,
    id.species =  idsp,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = file.path(tdir, "B_625_625_5spc_repl1.rds")
  )
  
  expect_equal(checkDataObject(testdat), "Data object looks ok.")
  
  tmp <- readRDS(normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"))
  expect_is(tmp, "phyloCompData")
  expect_is(count.matrix(tmp), "matrix")
  expect_equal(nrow(count.matrix(tmp)), 489) # some are filtered
  expect_equal(ncol(count.matrix(tmp)), 8)
  expect_equal(nrow(length.matrix(tmp)), 489) # some are filtered
  expect_equal(ncol(length.matrix(tmp)), 8)
  expect_is(sample.annotations(tmp), "data.frame")
  expect_equal(nrow(sample.annotations(tmp)), 8)
  expect_equal(ncol(sample.annotations(tmp)), 4)
  expect_is(variable.annotations(tmp), "data.frame")
  expect_equal(nrow(variable.annotations(tmp)), 489)
  expect_equal(ncol(variable.annotations(tmp)), 22)
  expect_is(info.parameters(tmp), "list")
  expect_equal(info.parameters(tmp)$n.diffexp, 50)
  expect_equal(info.parameters(tmp)$dataset, "B_625_625")
  expect_equal(info.parameters(tmp)$fraction.upregulated, 0.5)
  expect_equal(info.parameters(tmp)$between.group.diffdisp, FALSE)
  expect_equal(info.parameters(tmp)$filter.threshold.total, 1)
  expect_equal(info.parameters(tmp)$filter.threshold.mediancpm, 0)
  expect_equal(info.parameters(tmp)$fraction.non.overdispersed, 0)
  expect_equal(info.parameters(tmp)$random.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$random.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.high.prob, 0)
  expect_equal(info.parameters(tmp)$single.outlier.low.prob, 0)
  expect_equal(info.parameters(tmp)$effect.size, 1.5)
  expect_equal(info.parameters(tmp)$samples.per.cond, 4)
  expect_equal(info.parameters(tmp)$repl.id, 1)
  expect_equal(info.parameters(tmp)$seqdepth, 1e5)
  expect_equal(info.parameters(tmp)$minfact, 0.7)
  expect_equal(info.parameters(tmp)$maxfact, 1.4)
  expect_equal(info.parameters(tmp)$nEff, 0.3636364, tolerance = 1e-7)
  expect_equal(info.parameters(tmp)$nEffRatio, info.parameters(tmp)$nEff) # n/4 = 1
  expect_equal(filtering(tmp), "total count >= 1 ;  median cpm >= 0")
  expect_is(analysis.date(tmp), "character")
  expect_equal(analysis.date(tmp), "")
  expect_is(package.version(tmp), "character")
  expect_equal(package.version(tmp), "")
  expect_is(method.names(tmp), "list")
  expect_equal(method.names(tmp), list())
  
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    expect_message(
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "DESeq2.length",
      Rmdfunction = "DESeq2.length.createRmd",
      output.directory = tdir, fit.type = "parametric",
      test = "Wald", beta.prior = TRUE,
      independent.filtering = TRUE, cooks.cutoff = TRUE,
      impute.outliers = TRUE,
      extra.design.covariates = NULL,
      nas_as_ones = TRUE
    ),
    "all NAs in adjusted p values are replaced by 1"
    )
  }
  # limma is in Imports
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "lengthNorm.limma",
    Rmdfunction = "lengthNorm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM",
    extra.design.covariates = NULL,
    length.normalization = "TPM",
    data.transformation = "log2",
    block.factor = NULL
  )
  if (requireNamespace("statmod", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
      result.extent = "lengthNorm.limma.cor",
      Rmdfunction = "lengthNorm.limma.createRmd",
      output.directory = tdir, norm.method = "TMM",
      extra.design.covariates = NULL,
      length.normalization = "TPM",
      data.transformation = "log2",
      block.factor = "id.species"
    )
  }
  # phylolm is in Imports
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "phylolm_cpm",
    Rmdfunction = "phylolm.createRmd",
    output.directory = tdir, norm.method = "TMM",
    model = "OUfixedRoot", measurement_error = TRUE,
    extra.design.covariates = NULL,
    length.normalization = "none",
    data.transformation = "sqrt"
  )
  # phylolm is in Imports
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "phylolm",
    Rmdfunction = "phylolm.createRmd",
    output.directory = tdir, norm.method = "TMM",
    model = "BM", measurement_error = TRUE,
    extra.design.covariates = NULL,
    length.normalization = "TPM",
    data.transformation = "log2"
  )
  
  # with extra factor
  tmp <- readRDS(normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"))
  sample.annotations(tmp)$test_reg <- rnorm(nrow(sample.annotations(tmp)))
  sample.annotations(tmp)$test_fac <- factor(sample(c(0, 1)))
  saveRDS(tmp, normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"))
  
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    runDiffExp(
      data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"),
      result.extent = "DESeq2.length.factor",
      Rmdfunction = "DESeq2.length.createRmd",
      output.directory = tdir, fit.type = "parametric",
      test = "Wald", beta.prior = TRUE,
      independent.filtering = TRUE, cooks.cutoff = TRUE,
      impute.outliers = FALSE,
      extra.design.covariates = c("test_reg", "test_fac")
    )
  }
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "lengthNorm.limma.factor",
    Rmdfunction = "lengthNorm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM",
    extra.design.covariates = c("test_reg", "test_fac"),
    trend = TRUE,
    length.normalization = "TPM",
    data.transformation = "log2",
    block.factor = NULL
  )
  runDiffExp(
    data.file = normalizePath(file.path(tdir, "B_625_625_5spc_repl1.rds"), winslash = "/"), 
    result.extent = "phylolm.factor",
    Rmdfunction = "phylolm.createRmd",
    output.directory = tdir, norm.method = "TMM",
    model = "BM", measurement_error = FALSE,
    extra.design.covariates = c("test_reg", "test_fac"),
    length.normalization = "TPM",
    data.transformation = "log2"
  )
  
  methods <- c("DESeq2.length", "DESeq2.length.factor", "lengthNorm.limma", "lengthNorm.limma.cor", "lengthNorm.limma.factor", "phylolm.factor", "phylolm_cpm", "phylolm")
  
  ## Test show() method
  m <- "lengthNorm.limma"
  tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")),
                               winslash = "/"))
  show(tmp)
  count.matrix(tmp) <- count.matrix(tmp)[, 1:4]
  show(tmp)
  
  for (m in methods) {
    if (!(m %in% c("DESeq2.length", "DESeq2.length.factor")) || requireNamespace("DESeq2", quietly = TRUE)) {
      if (!(m %in% c("lengthNorm.limma.cor")) || requireNamespace("statmod", quietly = TRUE)) {
        tmp <- readRDS(normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), winslash = "/"))
        
        expect_is(tmp, "phyloCompData")
        expect_is(result.table(tmp), "data.frame")
        expect_equal(nrow(result.table(tmp)), 489)
        expect_is(code(tmp), "character")
        expect_is(analysis.date(tmp), "character")
        expect_is(compcodeR:::package.version(tmp), "character")
        expect_is(method.names(tmp), "list")
        expect_named(method.names(tmp), c("short.name", "full.name"))
        
        tmp2 <- tmp; result.table(tmp2) <- result.table(tmp2)[1:10, ]; expect_equal(check_compData_results(tmp2), "result.table must have the same number of rows as count.matrix.")
        tmp2 <- tmp; result.table(tmp2) <- data.frame(); expect_equal(check_compData_results(tmp2), "Object must contain a data frame named 'result.table'.")
        tmp2 <- tmp; result.table(tmp2)$score <- NULL; expect_equal(check_compData_results(tmp2), "result.table must contain a column named 'score'.")
      }
    }
  }
  
  for (m in methods) {
    if (!(m %in% c("DESeq2.length", "DESeq2.length.factor")) || requireNamespace("DESeq2", quietly = TRUE)) {
      if (!(m %in% c("lengthNorm.limma.cor")) || requireNamespace("statmod", quietly = TRUE)) {
        generateCodeHTMLs(
          normalizePath(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), 
                        winslash = "/"), normalizePath(tdir)
        )
        expect_true(file.exists(normalizePath(file.path(
          tdir, paste0("B_625_625_5spc_repl1_", 
                       m, "_code.html")), winslash = "/")))
      }
    }
  }
  
  ## Comparison report
  file.table <- data.frame(input.files = normalizePath(file.path(
    tdir, paste0("B_625_625_5spc_repl1_", 
                 methods[!(methods %in% c("DESeq2.length", "DESeq2.length.factor", "lengthNorm.limma.cor"))],
                 ".rds")), winslash = "/"))
  parameters <- NULL
  
  comp <- runComparison(file.table = file.table, output.directory = tdir,
                        parameters = parameters, save.result.table = TRUE, knit.results = FALSE)
  ff <- list.files(path = tdir, pattern = "compcodeR_result_table_.*.rds", full.names = TRUE)
  resTable <- readRDS(ff[length(ff)])
  expect_equal(resTable$fp + resTable$tp + resTable$fn + resTable$tn, rep(489, nrow(resTable)))
  expect_equal(ncol(resTable), 14)
  expect_equal(nrow(resTable), length(methods) - 3)
  
  parameters <- list()
  par2 <- parameters; par2$incl.dataset <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with datasets")
  par2 <- parameters; par2$incl.nbr.samples <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with nbr.samples")
  par2 <- parameters; par2$incl.replicates <- 10
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with replicates")
  par2 <- parameters; par2$incl.de.methods <- "missing"
  expect_error(runComparison(file.table = file.table, output.directory = tdir,
                             parameters = par2),
               "No methods left to compare after matching with DE methods")
})
