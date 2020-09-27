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
  tmp
  
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
  expect_equal(check_compData(tmp@count.matrix), "This is not an S4 object.")
  tmp2 <- tmp; tmp2@count.matrix <- as.matrix(numeric(0)); expect_equal(check_compData(tmp2), "Object must contain a non-empty count matrix.")
  tmp2 <- tmp; tmp2@sample.annotations <- as.data.frame(NULL); expect_equal(check_compData(tmp2), "Object must contain a non-empty sample annotation data frame.")
  tmp2 <- tmp; tmp2@sample.annotations <- data.frame(condition = numeric(0)); expect_equal(check_compData(tmp2), "The sample.annotations must contain a column named condition.")
  tmp2 <- tmp; tmp2@info.parameters <- list(); expect_equal(check_compData(tmp2), "Object must contain a non-empty list called info.parameters.")
  tmp2 <- tmp; tmp2@info.parameters <- list(tmp = 1); expect_equal(check_compData(tmp2), "info.parameters list must contain an entry named 'dataset'.")
  tmp2 <- tmp; tmp2@info.parameters <- list(dataset = ""); expect_equal(check_compData(tmp2), "info.parameters list must contain an entry named 'uID'.")
  tmp2 <- tmp; tmp2@count.matrix <- tmp2@count.matrix[1:10, ]; expect_equal(check_compData(tmp2), "count.matrix and variable.annotations do not contain the same number of rows.")
  tmp2 <- tmp; rownames(tmp2@count.matrix) <- paste0("r", rownames(tmp2@count.matrix)); expect_equal(check_compData(tmp2), "The rownames of count.matrix and variable.annotations are not the same.")
  tmp2 <- tmp; tmp2@count.matrix <- tmp2@count.matrix[, 1:2]; expect_equal(check_compData(tmp2), "The number of columns of count.matrix is different from the number of rows of sample.annotations.")
  tmp2 <- tmp; colnames(tmp2@count.matrix) <- paste0("r", colnames(tmp2@count.matrix)); expect_equal(check_compData(tmp2), "The colnames of count.matrix are different from the rownames of sample.annotations.")
  
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
  expect_equal(tmp@variable.annotations$truedispersions.S1, 
               tmp@variable.annotations$truedispersions.S2)
  expect_equal(tmp@variable.annotations$truemeans.S1,
               tmp@variable.annotations$truemeans.S2)
  expect_equal(sum(tmp@variable.annotations$n.random.outliers.up.S1 + 
                     tmp@variable.annotations$n.random.outliers.up.S2 + 
                     tmp@variable.annotations$n.random.outliers.down.S1 + 
                     tmp@variable.annotations$n.random.outliers.down.S2 + 
                     tmp@variable.annotations$n.single.outliers.up.S1 + 
                     tmp@variable.annotations$n.single.outliers.up.S2 + 
                     tmp@variable.annotations$n.single.outliers.down.S1 + 
                     tmp@variable.annotations$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(tmp@variable.annotations$truelog2foldchanges)), 0)
  expect_equal(sum(tmp@variable.annotations$upregulation + 
                     tmp@variable.annotations$downregulation + 
                     tmp@variable.annotations$differential.expression), 0)
               
  ## Specify effect sizes individually
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    effect.size = c(exp(abs(rnorm(5))), exp(-abs(rnorm(5))), rep(1, 40)),
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(tmp@variable.annotations$truedispersions.S1, 
               tmp@variable.annotations$truedispersions.S2)
  expect_equal(tmp@variable.annotations$truemeans.S1[-(1:10)],
               tmp@variable.annotations$truemeans.S2[-(1:10)])
  expect_equal(sum(tmp@variable.annotations$n.random.outliers.up.S1 + 
                     tmp@variable.annotations$n.random.outliers.up.S2 + 
                     tmp@variable.annotations$n.random.outliers.down.S1 + 
                     tmp@variable.annotations$n.random.outliers.down.S2 + 
                     tmp@variable.annotations$n.single.outliers.up.S1 + 
                     tmp@variable.annotations$n.single.outliers.up.S2 + 
                     tmp@variable.annotations$n.single.outliers.down.S1 + 
                     tmp@variable.annotations$n.single.outliers.down.S2), 0)
  expect_equal(sum(abs(tmp@variable.annotations$truelog2foldchanges[-(1:10)])), 0)
  expect_equal(sign(tmp@variable.annotations$truelog2foldchanges),
               c(rep(1, 5), rep(-1, 5), rep(0, 40)))
  expect_equal(tmp@variable.annotations$upregulation, 
               c(rep(1, 5), rep(0, 45)))
  expect_equal(tmp@variable.annotations$downregulation, 
               c(rep(0, 5), rep(1, 5), rep(0, 40)))
  expect_equal(sum(tmp@variable.annotations$upregulation + 
                     tmp@variable.annotations$downregulation), 10)
  
  
  ## Different dispersions between groups
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    between.group.diffdisp = TRUE,
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(sign(abs(tmp@variable.annotations$truedispersions.S1 -  
                        tmp@variable.annotations$truedispersions.S2)),
               rep(1, 50))
  expect_equal(tmp@variable.annotations$upregulation, 
               c(rep(1, 10), rep(0, 40)))
  expect_equal(tmp@variable.annotations$downregulation, 
               rep(0, 50))
  expect_equal(sum(tmp@variable.annotations$upregulation + 
                     tmp@variable.annotations$downregulation), 10)
  
  ## Not overdispersed
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    between.group.diffdisp = FALSE, 
    fraction.non.overdispersed = 0.5,
    output.file = NULL
  )
  expect_is(tmp, "compData")
  expect_equal(tmp@variable.annotations$truedispersions.S1,  
               tmp@variable.annotations$truedispersions.S2)
  expect_equal(any(tmp@variable.annotations$truedispersions.S1 == 0), TRUE)
  expect_equal(tmp@info.parameters$fraction.non.overdispersed, 0.5)

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
  expect_equal(any(tmp@variable.annotations$n.random.outliers.up.S1 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.random.outliers.up.S2 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.random.outliers.down.S1 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.random.outliers.down.S2 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.single.outliers.up.S1 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.single.outliers.up.S2 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.single.outliers.down.S1 > 0), TRUE)
  expect_equal(any(tmp@variable.annotations$n.single.outliers.down.S2 > 0), TRUE)
  
  ## Summary report
  expect_error(summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp.rds")),
               "output.file must be an .html file.")
  summarizeSyntheticDataSet(tmp, file.path(tdir, "tmp_summaryrep.html"))
  expect_equal(file.exists(file.path(tdir, "tmp_summaryrep.html")), TRUE)
})

test_that("help functions work", {
  listcreateRmd()
  
  tmp <- generateSyntheticData(
    dataset = "B_625_625", n.vars = 50, 
    samples.per.cond = 5, n.diffexp = 10,
    output.file = NULL
  )
  
  expect_error(
    computeMval(tmp@count.matrix, c(3, tmp@sample.annotations$condition[-1])),
    "Must have exactly two groups to calculate M-value"
  )
  expect_error(
    computeAval(tmp@count.matrix, c(3, tmp@sample.annotations$condition[-1])),
    "Must have exactly two groups to calculate A-value"
  )
  
  mval <- computeMval(tmp@count.matrix, tmp@sample.annotations$condition)
  aval <- computeAval(tmp@count.matrix, tmp@sample.annotations$condition)
  
  expect_is(mval, "numeric")
  expect_is(aval, "numeric")
  expect_equal(length(mval), nrow(tmp@count.matrix))
  expect_equal(length(aval), nrow(tmp@count.matrix))
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
  expect_is(tmp@count.matrix, "matrix")
  expect_equal(nrow(tmp@count.matrix), 500)
  expect_equal(ncol(tmp@count.matrix), 10)
  expect_is(tmp@sample.annotations, "data.frame")
  expect_equal(nrow(tmp@sample.annotations), 10)
  expect_equal(ncol(tmp@sample.annotations), 2)
  expect_is(tmp@variable.annotations, "data.frame")
  expect_equal(nrow(tmp@variable.annotations), 500)
  expect_equal(ncol(tmp@variable.annotations), 18)
  expect_is(tmp@info.parameters, "list")
  expect_equal(tmp@info.parameters$n.diffexp, 50)
  expect_equal(tmp@info.parameters$dataset, "B_625_625")
  expect_equal(tmp@info.parameters$fraction.upregulated, 0.5)
  expect_equal(tmp@info.parameters$between.group.diffdisp, FALSE)
  expect_equal(tmp@info.parameters$filter.threshold.total, 1)
  expect_equal(tmp@info.parameters$filter.threshold.mediancpm, 0)
  expect_equal(tmp@info.parameters$fraction.non.overdispersed, 0)
  expect_equal(tmp@info.parameters$random.outlier.high.prob, 0)
  expect_equal(tmp@info.parameters$random.outlier.low.prob, 0)
  expect_equal(tmp@info.parameters$single.outlier.high.prob, 0)
  expect_equal(tmp@info.parameters$single.outlier.low.prob, 0)
  expect_equal(tmp@info.parameters$effect.size, 1.5)
  expect_equal(tmp@info.parameters$samples.per.cond, 5)
  expect_equal(tmp@info.parameters$repl.id, 1)
  expect_equal(tmp@info.parameters$seqdepth, 1e5)
  expect_equal(tmp@info.parameters$minfact, 0.7)
  expect_equal(tmp@info.parameters$maxfact, 1.4)
  expect_equal(tmp@filtering, "total count >= 1 ;  median cpm >= 0")
  expect_is(tmp@analysis.date, "character")
  expect_equal(tmp@analysis.date, "")
  expect_is(tmp@package.version, "character")
  expect_equal(tmp@package.version, "")
  expect_is(tmp@method.names, "list")
  expect_equal(tmp@method.names, list())

  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "baySeq",
    Rmdfunction = "baySeq.createRmd",
    output.directory = tdir, norm.method = "edgeR",
    equaldisp = TRUE
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "DESeq2",
    Rmdfunction = "DESeq2.createRmd",
    output.directory = tdir, fit.type = "parametric",
    test = "Wald", beta.prior = TRUE,
    independent.filtering = TRUE, cooks.cutoff = TRUE,
    impute.outliers = TRUE
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "DSS",
    Rmdfunction = "DSS.createRmd",
    output.directory = tdir, norm.method = "quantile",
    disp.trend = FALSE
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "EBSeq",
    Rmdfunction = "EBSeq.createRmd",
    output.directory = tdir, norm.method = "median"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "edgeR.exact",
    Rmdfunction = "edgeR.exact.createRmd",
    output.directory = tdir, norm.method = "TMM",
    trend.method = "movingave", disp.type = "tagwise"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "edgeR.GLM",
    Rmdfunction = "edgeR.GLM.createRmd",
    output.directory = tdir, norm.method = "TMM",
    disp.type = "tagwise", disp.method = "CoxReid",
    trended = TRUE
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "logcpm.limma",
    Rmdfunction = "logcpm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "NBPSeq",
    Rmdfunction = "NBPSeq.createRmd",
    output.directory = tdir, norm.method = "TMM",
    disp.method = "NBP"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "NOISeq",
    Rmdfunction = "NOISeq.prenorm.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "sqrtcpm.limma",
    Rmdfunction = "sqrtcpm.limma.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "TCC",
    Rmdfunction = "TCC.createRmd",
    output.directory = tdir, norm.method = "tmm",
    test.method = "edger", iteration = 3,
    normFDR = 0.1, floorPDEG = 0.05
  )
  runDiffExp(
    data.file = file.path(tdir, "B_625_625_5spc_repl1.rds"), 
    result.extent = "ttest",
    Rmdfunction = "ttest.createRmd",
    output.directory = tdir, norm.method = "TMM"
  )
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
               "voom.ttest")

  for (m in methods) {
    tmp <- readRDS(file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")))
    
    expect_is(tmp, "compData")
    expect_is(tmp@result.table, "data.frame")
    expect_equal(nrow(tmp@result.table), 500)
    expect_is(tmp@code, "character")
    expect_is(tmp@analysis.date, "character")
    expect_is(tmp@package.version, "character")
    expect_is(tmp@method.names, "list")
    expect_named(tmp@method.names, c("short.name", "full.name"))
    
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

})
