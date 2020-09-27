test_that("runComparison works", {
  tdir <- tempdir()
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
  }
  
  for (m in methods) {
    generateCodeHTMLs(
      file.path(tdir, paste0("B_625_625_5spc_repl1_", m, ".rds")), tdir
    )
    expect_true(file.exists(file.path(tdir, paste0("B_625_625_5spc_repl1_", 
                                                   m, "_code.html"))))
  }

})
