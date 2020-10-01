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

test_that("generateSyntheticData works", {
  tdir <- tempdir()
  
  ## No DEGs
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
  expect_equal(file.exists(file.path(tdir, "tmp_summaryrep.html")), TRUE)
})

test_that("help functions work", {
  listcreateRmd()
  
  expect_error(checkRange("hello", "name", 0, 1), "Illegal value")
  expect_equal(checkRange(-1, "name", 0, 1), 0)
  expect_equal(checkRange(2, "name", 0, 1), 1)
  expect_equal(checkRange("-1", "name", 0, 1), 0)
  
  expect_equal(shorten.method.names(c("AUC", "ROC, all replicates")),
               c("auc", "rocall"))
  
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
  
  tmp <- readRDS(file.path(tdir, "B_625_625_5spc_repl1.rds"))
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

  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "baySeq",
  #   Rmdfunction = "baySeq.createRmd",
  #   output.directory = tdir, norm.method = "edgeR",
  #   equaldisp = TRUE
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "DESeq2",
  #   Rmdfunction = "DESeq2.createRmd",
  #   output.directory = tdir, fit.type = "parametric",
  #   test = "Wald", beta.prior = TRUE,
  #   independent.filtering = TRUE, cooks.cutoff = TRUE,
  #   impute.outliers = TRUE
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "DSS",
  #   Rmdfunction = "DSS.createRmd",
  #   output.directory = tdir, norm.method = "quantile",
  #   disp.trend = FALSE
  # )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"),
    result.extent = "EBSeq",
    Rmdfunction = "EBSeq.createRmd",
    output.directory = tdir, norm.method = "median"
  )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "edgeR.exact",
  #   Rmdfunction = "edgeR.exact.createRmd",
  #   output.directory = tdir, norm.method = "TMM",
  #   trend.method = "movingave", disp.type = "tagwise"
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "edgeR.GLM",
  #   Rmdfunction = "edgeR.GLM.createRmd",
  #   output.directory = tdir, norm.method = "TMM",
  #   disp.type = "tagwise", disp.method = "CoxReid",
  #   trended = TRUE
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "logcpm.limma",
  #   Rmdfunction = "logcpm.limma.createRmd",
  #   output.directory = tdir, norm.method = "TMM"
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "NBPSeq",
  #   Rmdfunction = "NBPSeq.createRmd",
  #   output.directory = tdir, norm.method = "TMM",
  #   disp.method = "NBP"
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "NOISeq",
  #   Rmdfunction = "NOISeq.prenorm.createRmd",
  #   output.directory = tdir, norm.method = "TMM"
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "sqrtcpm.limma",
  #   Rmdfunction = "sqrtcpm.limma.createRmd",
  #   output.directory = tdir, norm.method = "TMM"
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "TCC",
  #   Rmdfunction = "TCC.createRmd",
  #   output.directory = tdir, norm.method = "tmm",
  #   test.method = "edger", iteration = 3,
  #   normFDR = 0.1, floorPDEG = 0.05
  # )
  # runDiffExp(
  #   data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
  #   result.extent = "ttest",
  #   Rmdfunction = "ttest.createRmd",
  #   output.directory = tdir, norm.method = "TMM"
  # )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "voom.limma",
    Rmdfunction = "voom.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "voom.ttest",
    Rmdfunction = "voom.ttest.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )

  methods <- c("baySeq", "DESeq2", "DSS", "EBSeq", "edgeR.exact",
               "edgeR.GLM", "logcpm.limma", "NBPSeq", "NOISeq", 
               "sqrtcpm.limma", "TCC", "ttest", "voom.limma",
               "voom.ttest")[c(4, 13, 14)]

  ## Test show() method
  m <- "EBSeq"
  tmp <- readRDS(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")))
  show(tmp)
  count.matrix(tmp) <- count.matrix(tmp)[, 1:4]
  show(tmp)
  
  for (m in methods) {
    tmp <- readRDS(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")))
    
    expect_is(tmp, "compData")
    expect_is(result.table(tmp), "data.frame")
    expect_equal(nrow(result.table(tmp)), 500)
    expect_is(code(tmp), "character")
    expect_is(analysis.date(tmp), "character")
    expect_is(package.version(tmp), "character")
    expect_is(method.names(tmp), "list")
    expect_named(method.names(tmp), c("short.name", "full.name"))
    
    tmp2 <- tmp; result.table(tmp2) <- result.table(tmp2)[1:10, ]; expect_equal(check_compData_results(tmp2), "result.table must have the same number of rows as count.matrix.")
    tmp2 <- tmp; result.table(tmp2) <- data.frame(); expect_equal(check_compData_results(tmp2), "Object must contain a data frame named 'result.table'.")
    tmp2 <- tmp; result.table(tmp2)$score <- NULL; expect_equal(check_compData_results(tmp2), "result.table must contain a column named 'score'.")
  }
  
  for (m in methods) {
    generateCodeHTMLs(
      file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), tdir
    )
    expect_true(file.exists(file.path(tdir, paste0("B_625_625_5spc_repl1_", 
                                                   m, "_code.html"))))
  }

  ## Comparison report
  file.table <- data.frame(input.files = file.path(
    tdir, paste0("B_625_625_5spc_repl1_", 
                 c("voom.limma", "voom.ttest", "EBSeq"), ".rds")))
  parameters <- NULL
  comp <- runComparison(file.table = file.table, output.directory = tdir,
                        parameters = parameters)
  
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
