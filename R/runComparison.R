# Generate a .Rmd file containing code to perform a comparison of differential
# expression methods
# 
# A function to generate code that can be run to perform a comparison of the
# performance of several differential expression methods. The code is written to
# a .Rmd file. This function is generally not called by the user, the main
# interface for comparing differential expression methods is the
# \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param output.file The path to the file where the code will be written.
# @author Charlotte Soneson
createResultsRmdFile <- function(setup.parameters.file, output.file) {
  ## Check that the output file ends with .Rmd
  if (!(substr(output.file, nchar(output.file) - 3, 
               nchar(output.file)) == ".Rmd")) {
    output.file <- sub(strsplit(output.file, 
                                "\\.")[[1]][length(strsplit(output.file, 
                                                            "\\.")[[1]])], 
                       "Rmd", 
                       output.file)
  }
  
  ## Open the output file
  resultfile <- file(output.file, open = "w")
  
  ## Load the setup parameters
  setup.parameters <- readRDS(setup.parameters.file)
  setup.parameters$incl.nbr.samples <- 
    as.numeric(setup.parameters$incl.nbr.samples)
  setup.parameters$incl.replicates <- 
    as.numeric(setup.parameters$incl.replicates)
  
  if (all(!is.na(setup.parameters$incl.nbr.samples))) {
    setup.parameters$incl.nbr.samples <- 
      sort(setup.parameters$incl.nbr.samples)
  }
  if (all(!is.na(setup.parameters$incl.replicates))) {
    setup.parameters$incl.replicates <- 
      sort(setup.parameters$incl.replicates)
  }
  
  ## Write the text that goes in the top of of the report (detailing which
  ## parameter values were used etc.)
  writeLines(c("```{r setup, echo = FALSE}",
               "options(width = 80)",
               "```"), resultfile)
  writeLines("# Comparison report, differential expression of RNAseq data", resultfile)
  writeLines(c(paste("Created by the compcodeR package, version", 
                     packageVersion("compcodeR")),
               paste("Date:", date())
               ##               paste("Parameter file:", setup.parameters.file), 
               ##               "```{r loadparam, echo = FALSE}",
               ##               paste("setup.parameters = readRDS('", setup.parameters.file, "')", sep = ""),
               ##               "```",
               ##paste("Data set:", setup.parameters$incl.dataset)
  ), resultfile, sep = "\n\n")
  writeLines("Data set:", resultfile)
  writeLines(c("```{r dataset, echo = FALSE}",
               "cat(setup.parameters$incl.dataset)",
               "```"), resultfile, sep = '\n')
  
  if (any(!is.na(setup.parameters$incl.nbr.samples))) {
    writeLines("Number of samples per condition:", resultfile)
    writeLines(c("```{r nsamples, echo = FALSE}",
                 "cat(sort(setup.parameters$incl.nbr.samples), sep = ', ')",
                 "```"), resultfile, sep = '\n')
  }
  if (any(!is.na(setup.parameters$incl.replicates))) {
    writeLines("Included replicates (for repeated simulated data sets):", 
               resultfile)
    writeLines(c("```{r replicates, echo = FALSE}", 
                 "cat(sort(setup.parameters$incl.replicates), sep = ', ')", 
                 "```"), resultfile, sep = '\n')
  }
  writeLines("Differential expression methods included in the comparison:", 
             resultfile)
  writeLines(c("```{r demethods, echo = FALSE}",
               "cat(setup.parameters$incl.de.methods, sep = '\n')",
               "```"), resultfile)
  
  writeLines("Parameter values:", resultfile)
  writeLines(c("```{r parametervalues, echo = FALSE}",
               "print(data.frame(value = unlist(setup.parameters[c('fdr.threshold', 'tpr.threshold', 'typeI.threshold', 'mcc.threshold', 'ma.threshold', 'fdc.maxvar', 'overlap.threshold', 'fracsign.threshold', 'nbrtpfp.threshold', 'signal.measure')])))",
               "```"), resultfile)
  
  writeLines("---", resultfile, sep = "\n\n")
  
  ## Set some parameters for the plotting
  nbr.cols <- 3#min(3, length(setup.parameters$incl.replicates) + 1)
  nbr.rows <- ceiling((length(setup.parameters$incl.replicates) + 1)/nbr.cols)
  nbr.cols2 <- 3#min(3, length(setup.parameters$incl.de.methods))
  nbr.rows2 <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols2)
  kk1 <- length(setup.parameters$incl.de.methods) * 
    length(setup.parameters$incl.nbr.samples)
  
  writeLines("<a name='contents'></a>", resultfile)
  writeLines("## Contents", resultfile)
  if ('rocall' %in% setup.parameters$comparisons) 
    writeLines("- [ROC curves, all replicates](#rocall)\n", resultfile)
  if ('rocone' %in% setup.parameters$comparisons)
    writeLines("- [ROC curves, single replicate](#rocone)\n", resultfile)
  if ('auc' %in% setup.parameters$comparisons)
    writeLines("- [Area under the ROC curve](#auc)\n", resultfile)
  if ('mcc' %in% setup.parameters$comparisons)
    writeLines("- [Matthew's correlation coefficient](#mcc)\n", resultfile)
  if ('fdcurvesall' %in% setup.parameters$comparisons)
    writeLines("- [False discovery curves, all replicates](#fdcurvesall)\n", 
               resultfile)
  if ('fdcurvesone' %in% setup.parameters$comparisons)
    writeLines("- [False discovery curves, single replicate](#fdcurvesone)\n", 
               resultfile)
  if ('maplot' %in% setup.parameters$comparisons)
    writeLines("- [MA plots, single replicate](#maplot)\n", resultfile)
  if ('scorevsexpr' %in% setup.parameters$comparisons)
    writeLines("- [Gene score vs average expression level, single replicate](#scorevsexpr)\n", resultfile)
  if ('scorevssignal' %in% setup.parameters$comparisons)
    writeLines("- [Gene score vs 'signal' for genes expressed in only one condition, all replicates](#scorevssignal)\n", resultfile)
  if ('scorevsoutlier' %in% setup.parameters$comparisons)
    writeLines("- [Score distribution vs number of outliers, single replicate](#scorevsoutlier)\n", resultfile)
  if ('fracsign' %in% setup.parameters$comparisons)
    writeLines("- [Fraction significant genes](#fracsign)\n", resultfile)
  if ('nbrsign' %in% setup.parameters$comparisons)
    writeLines("- [Number of significant genes](#nbrsign)\n", resultfile)
  if ('nbrtpfp' %in% setup.parameters$comparisons)
    writeLines("- [Number of TP, FP, TN, FN](#nbrtpfp)\n", resultfile)
  if ('typeIerror' %in% setup.parameters$comparisons)
    writeLines("- [Type I error](#typeIerror)\n", resultfile)
  if ('fdr' %in% setup.parameters$comparisons)
    writeLines("- [False discovery rate](#fdr)\n", resultfile)
  if ('fdrvsexpr' %in% setup.parameters$comparisons)
    writeLines("- [False discovery rate vs average expression level](#fdrvsexpr)\n", resultfile)
  if ('tpr' %in% setup.parameters$comparisons)
    writeLines("- [True positive rate](#tpr)\n", resultfile)
  if ('overlap' %in% setup.parameters$comparisons)
    writeLines("- [Overlap between sets of differentially expressed genes, single replicate](#overlap)\n", resultfile)
  if ('sorensen' %in% setup.parameters$comparisons)
    writeLines("- [Sorensen index between sets of differentially expressed genes, single replicate](#sorensen)\n", resultfile)
  if ('correlation' %in% setup.parameters$comparisons)
    writeLines("- [Spearman correlation between scores, single replicate](#correlation)\n", resultfile)
  writeLines("---", resultfile)
  
  ## Include the actual results
  if ('rocall' %in% setup.parameters$comparisons) {
    writeLines("<a name='rocall'></a>", resultfile)
    writeLines("## ROC (receiver operating characteristic) curves [(Contents)](#contents)", resultfile)
    writeLines("A receiver operating characteristic (ROC) curve is a way to summarize the ability of a test or ranking procedure to rank truly positive (i.e., truly differentially expressed) genes ahead of truly negative (i.e., truly non-differentially expressed). To create the ROC curve for a given differential expression method, the genes are ranked in decreasing order by the score, which is assigned to them during the differential expression analysis and quantifies the degree of statistical significance or association with the predictor (the condition). For a given threshold, all genes with scores above the threshold are classified as 'positive' and all genes with scores below the threshold are classified as 'negative'. Comparing these assignments to the true differential expression status, a true positive rate and a false positive rate can be computed and marked in a plot. As the threshold is changed, these pairs of values trace out the ROC curve. A good test procedure gives a ROC curve which passes close to the upper left corner of the plot, while a poor test corresponds to a ROC curve close to the diagonal. \n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i],
                         " samples/condition [(Contents)](#contents)"), resultfile)
      }
      writeLines(c(paste("```{r rocall-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols, ", fig.height = ",
                         6*nbr.rows, ", fig.align = 'center'}", sep = ""),
                   paste("makeROCcurves(setup.parameters, sel.nbrsamples = ",
                   setup.parameters$incl.nbr.samples[i], ", sel.repl =  c(", paste(intersect(setup.parameters$incl.replicates, setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)]), collapse = ", "), "))", sep = ""),
                   "```", "---"), resultfile)
    }
  }
  
  if ('rocone' %in% setup.parameters$comparisons) {
    writeLines("<a name='rocone'></a>", resultfile)
    writeLines("## ROC curves, single replicate [(Contents)](#contents)", resultfile)
    writeLines("A receiver operating characteristic (ROC) curve is a way to summarize the ability of a test or ranking procedure to rank truly positive (i.e., truly differentially expressed) genes ahead of truly negative (i.e., truly non-differentially expressed). To create the ROC curve for a given differential expression method, the genes are ranked in decreasing order by the score, which is assigned to them during the differential expression analysis and quantifies the degree of statistical significance or association with the predictor (the condition). For a given threshold, all genes with scores above the threshold are classified as 'positive' and all genes with scores below the threshold are classified as 'negative'. Comparing these assignments to the true differential expression status, a true positive rate and a false positive rate can be computed and marked in a plot. As the threshold is changed, these pairs of values trace out the ROC curve. A good test procedure gives a ROC curve which passes close to the upper left corner of the plot, while a poor test corresponds to a ROC curve close to the diagonal. \n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i],
                         " samples/condition [(Contents)](#contents)"), resultfile)
      }
      writeLines(c(paste("```{r rocone-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 18, fig.height = 6, fig.align = 'left'}", sep = ""),
                   paste("makeROCcurves(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl = ",
                         mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")", sep = ""),
                   "```", "---"), resultfile)
    }
  }
  
  if (any(c('auc', 'fdr', 'tpr', 'typeIerror', 'fracsign', 'nbrtpfp', 'nbrsign', 'mcc') %in% setup.parameters$comparisons)) {
    writeLines(c(paste("```{r tabl, eval = TRUE, include = TRUE}"),
                 "res.table = createResultTable(setup.parameters)"), resultfile)
    if (any(!is.na(setup.parameters$incl.nbr.samples))) {
      writeLines("res.table = padResultTable(res.table)", resultfile)
    }
    writeLines(c("```", "---"), resultfile)
  }
  
  if ('auc' %in% setup.parameters$comparisons) {
    writeLines("<a name='auc'></a>", resultfile)
    writeLines("## AUC [(Contents)](#contents)", resultfile)
    writeLines("A receiver operating characteristic (ROC) curve is a way to summarize the ability of a test or ranking procedure to rank truly positive (i.e., truly differentially expressed) genes ahead of truly negative (i.e., truly non-differentially expressed). To create the ROC curve for a given differential expression method, the genes are ranked in decreasing order by the score, which is assigned to them during the differential expression analysis and quantifies the degree of statistical significance or association with the predictor (the condition). For a given threshold, all genes with scores above the threshold are classified as 'positive' and all genes with scores below the threshold are classified as 'negative'. Comparing these assignments to the true differential expression status, a true positive rate and a false positive rate can be computed and marked in a plot. As the threshold is changed, these pairs of values trace out the ROC curve. A good test procedure gives a ROC curve which passes close to the upper left corner of the plot, while a poor test corresponds to a ROC curve close to the diagonal. The area under the ROC curve (AUC) summarizes the performance of the ranking procedure. A good method gives an AUC close to 1, while a poor method gives an AUC closer to 0.5. Each boxplot below summarizes the AUCs across all data set replicates included in the comparison. \n", resultfile)
    writeLines(c(paste("```{r auc, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height = ",
                       length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'auc')",
                 "```", "---"), resultfile)
  }
  
  if ('mcc' %in% setup.parameters$comparisons) {
    writeLines("<a name='mcc'></a>", resultfile)
    writeLines("## Matthew's correlation coefficient [(Contents)](#contents)", resultfile)
    writeLines(paste0("Matthew's correlation coefficient is a measure summarizing the performance of a binary classifier, by combining the number of observed true positives, true negatives, false positives and false negatives into a single value. The correlation coefficient takes values between -1 and 1, where 1 corresponds to perfect classification and -1 corresponds to complete disagreement between the classification and the true labels. A value of 0 corresponds to random assignment. The figures below indicate the MCC obtained at an adjusted p-value threshold of ", setup.parameters$mcc.threshold, ". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n"), resultfile)
    writeLines(c(paste("```{r mcc, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height = ",
                       length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'mcc')",
                 "```", "---"), resultfile)
  }
  
  if ('fdcurvesall' %in% setup.parameters$comparisons) {
    writeLines("<a name='fdcurvesall'></a>", resultfile)
    writeLines("## False discovery curves, all replicates [(Contents)](#contents)", resultfile)
    writeLines("A false discovery curve depicts the number of false positives encountered while stepping through a list of genes ranked by a score representing their statistical significance or the degree of association with a predictor. The truly differentially expressed genes are considered the 'true positives', and the truly non-differentially expressed genes the 'true negatives'. Hence, at a given position in the ranked list (shown on the x-axis) the value on the y-axis represents the number of truly non-differentially expressed genes that are ranked above that position. A good ranking method puts few true negatives among the top-ranked genes, and hence the false discovery curve rises slowly. A poor ranking method is recognized by a steeply increasing false discovery curve.\n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r fdcall-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols, ", fig.height = ", 6*nbr.rows, ", fig.align = 'left'}", sep = ""),
                   paste("makeFalseDiscoveryCurves(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl =  c(", paste(intersect(setup.parameters$incl.replicates, setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)]), collapse = ", "), "))", sep = ""),
                   "```", "---"), resultfile)
    }
  }
  
  if ('fdcurvesone' %in% setup.parameters$comparisons) {
    writeLines("<a name='fdcurvesone'></a>", resultfile)
    writeLines("## False discovery curves, single replicate [(Contents)](#contents)", resultfile)
    writeLines("A false discovery curve depicts the number of false positives encountered while stepping through a list of genes ranked by a score representing their statistical significance or the degree of association with a predictor. The truly differentially expressed genes are considered the 'true positives', and the truly non-differentially expressed genes the 'true negatives'. Hence, at a given position in the ranked list (shown on the x-axis) the value on the y-axis represents the number of truly non-differentially expressed genes that are ranked above that position. A good ranking method puts few true negatives among the top-ranked genes, and hence the false discovery curve rises slowly. A poor ranking method is recognized by a steeply increasing false discovery curve.\n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r fdcone-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 18, fig.height = 6, fig.align = 'left'}", sep = ""),
                   paste("makeFalseDiscoveryCurves(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl = ",
                         mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")"),
                   "```", "---"), resultfile)
    }
  }
  
  if ('maplot' %in% setup.parameters$comparisons) {
    writeLines("<a name='maplot'></a>", resultfile)
    writeLines("## MA plots, single replicate [(Contents)](#contents)", resultfile)
    writeLines(paste("An MA plot depicts the average expression level of the genes ('A', shown on the x-axis) and their log-fold change between two conditions ('M', shown on the y-axis). The genes called differentially expressed at an adjusted p-value threshold of", setup.parameters$ma.threshold, "are marked in color.\n"), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r maplot-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols2, ", fig.height = ", 6*nbr.rows2, ", fig.align = 'left'}", sep = ""),
                   paste("plotMASignificant(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl = ",
                         mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")"),
                   "```", "---"), resultfile)
    }
  }
  
  if ('scorevsexpr' %in% setup.parameters$comparisons) {
    writeLines("<a name='scorevsexpr'></a>", resultfile)
    writeLines("## Gene score vs average expression level, single replicate [(Contents)](#contents)", resultfile)
    writeLines(paste("In the figures below the gene score, which is computed in the differential expression analysis (and stored in the 'score' field of the result object), is plotted against the average expression level of the genes ('A', shown on the x-axis). A high value of the score signifies a 'more significant' gene. The colored line represents a loess fit to the data.\n"), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r scorevsexpr-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols2, ", fig.height = ", 6*nbr.rows2, ", fig.align = 'left'}", sep = ""),
                   paste("plotScoreVsExpr(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl = ",
                         mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")"),
                   "```", "---"), resultfile)
    }
  }
  
  if ('scorevssignal' %in% setup.parameters$comparisons) {
    writeLines("<a name='scorevssignal'></a>", resultfile)
    writeLines("## Gene score vs 'signal' for genes expressed in only one condition, all replicates [(Contents)](#contents)", resultfile)
    writeLines(paste("In the figures below the gene score, which is computed in the differential expression analysis (and stored in the 'score' field of the result object), is plotted against the 'signal', for genes that are expressed in only one of the two conditions. The signal (shown on the x-axis) is defined by computing the ", c("logarithm (base 2) of the normalized pseudo-counts for all samples in the condition where the gene is expressed, and averaging these values", "signal-to-noise ratio in the condition where the gene is expressed, i.e., first computing the logarithm (base 2) of the normalized pseudo-counts for all samples in this condition, and then computing the ratio between the mean and the standard deviation of these values")[match(setup.parameters$signal.measure, c("mean", "snr"))], ". A high value of the score (shown on the y-axis) signifies a gene that is more strongly associated with the predictor (the sample condition). We expect the methods to give higher scores to genes with stronger signal. The results from all dataset replicates included in the comparison are shown in the same plot, and the points are colored according to the data set replicate they originate from.\n"), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r scorevssignal-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols2, ", fig.height = ", 6*nbr.rows2, ", fig.align = 'left'}", sep = ""),
                   paste("plotSignalForZeroCounts(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl =  c(", paste(intersect(setup.parameters$incl.replicates, setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)]), collapse = ", "), "))", sep = ""),
                   "```", "---"), resultfile)
    }
  }
  
  if ('scorevsoutlier' %in% setup.parameters$comparisons) {
    writeLines("<a name='scorevsoutlier'></a>", resultfile)
    writeLines("## Score distribution vs number of outliers, single replicate [(Contents)](#contents)", resultfile)
    writeLines("The violin plots below show the distribution of the score assigned to the genes by the differential expression methods, as a function of the number of 'outlier counts' (that is, extremely high or low counts introduced artificially in the data and not generated by the underlying statistical distribution) for the genes. All types of outliers are summed. This allows an investigation of the sensitivity of a differential expression method to outlier counts (deviations from the underlying statistical model). A method that is sensitive to outliers shows a different score distribution for genes with outlier counts than for genes without outlier counts. When interpreting the figures below, be observant on the number of genes generating each distribution (indicated below the figure), since an empirical distribution based on only a few values many not be completely representative of the true distribution.\n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r scorevsoutlier-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols2, ", fig.height = ", 6*nbr.rows2, ", fig.align = 'left'}", sep = ""),
                   paste("plotScoreVsOutliers(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl = ",
                         mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")"),
                   "```", "---"), resultfile)
    }
  }
  
  
  if ('fracsign' %in% setup.parameters$comparisons) {
    writeLines("<a name='fracsign'></a>", resultfile)
    writeLines("## Fraction significant [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the fraction of genes in the data set that are called significant at an adjusted p-value threshold of ", setup.parameters$fracsign.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r fracsign, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'fracsign')",
                 "```", "---"), resultfile)
  }
  
  if ('nbrsign' %in% setup.parameters$comparisons) {
    writeLines("<a name='nbrsign'></a>", resultfile)
    writeLines("## Number of significant genes [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the number of genes in the data set that are called significant at an adjusted p-value threshold of ", setup.parameters$fracsign.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r nbrsign, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'nbrsign')",
                 "```", "---"), resultfile)
  }
  
  if ('nbrtpfp' %in% setup.parameters$comparisons) {
    writeLines("<a name='nbrtpfp'></a>", resultfile)
    writeLines("## Number of true positives [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the number of true positive genes found at an adjusted p-value threshold of ", setup.parameters$nbrtpfp.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r nbrtp, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'tp')",
                 "```", "---"), resultfile)
    
    writeLines("## Number of false positives [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the number of false positive genes found at an adjusted p-value threshold of ", setup.parameters$nbrtpfp.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r nbrfp, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'fp')",
                 "```", "---"), resultfile)
    
    writeLines("## Number of true negatives [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the number of true negative genes found at an adjusted p-value threshold of ", setup.parameters$nbrtpfp.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r nbrtn, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'tn')",
                 "```", "---"), resultfile)
    
    writeLines("## Number of false negatives [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below indicate the number of false negative genes found at an adjusted p-value threshold of ", setup.parameters$nbrtpfp.threshold,". Only differential expression methods returning corrected p-values or FDR estimates are included in the figure. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n",
                     sep = ''), resultfile)
    writeLines(c(paste("```{r nbrfn, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'fn')",
                 "```", "---"), resultfile)
  }
  
  if ('typeIerror' %in% setup.parameters$comparisons) {
    writeLines("<a name='typeIerror'></a>", resultfile)
    writeLines("## Type I error [(Contents)](#contents)", resultfile)
    writeLines(paste("The nominal p-value returned by a statistical test indicate the probability of obtaining a value of the test statistic which is at least as extreme as the one observed, given that the null hypothesis is true, e.g., that the gene is not differentially expressed between the compared conditions. Classifying a gene as statistically significantly differentially expressed although it is not truly differentially expressed is referred to as a type I error. For a good statistical test, we expect that the observed type I error rate at a given nominal p-value threshold (that is, the fraction of truly non-differentially expressed genes with a nominal p-value below this threshold) does not exceed the threshold. The figures below show the observed type I error rate at a nominal p-value threshold of ", setup.parameters$typeI.threshold, ". Only differential expression methods returning nominal p-values are included in the comparison. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n", sep = ''), resultfile)
    writeLines(c(paste("```{r typeIerror, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'typeIerror')",
                 "```", "---"), resultfile)
  }
  
  if ('fdr' %in% setup.parameters$comparisons) {
    writeLines("<a name='fdr'></a>", resultfile)
    writeLines("## FDR [(Contents)](#contents)", resultfile)
    writeLines(paste("The false discovery rate (FDR) indicates the fraction of truly non-differentially expressed genes that we expect to find among the genes that we consider to be differentially expressed. For high-dimensional problems, where many statistical tests are performed simultaneously (such as gene expression studies) it is more relevant to attempt to control the FDR than to control the gene-wise type I error rate, since it is almost certain that at least one gene will show a low nominal p-value even if the null hypothesis is true. To control the FDR, typically, the nominal p-values are adjusted for the large number of tests that are performed. The figures below indicate the observed rate of false discoveries (i.e., the fraction of truly non-differentially expressed genes among the genes that are considered significant) at an adjusted p-value threshold of ", setup.parameters$fdr.threshold,". Only methods returning corrected p-values or FDR estimates are included. Each boxplot summarizes the values obtained across all data set replicates that are included in the comparison. For a good method, the observed FDR should not be too high above the imposed adjusted p-value threshold (indicated by a dashed vertical line). If the observed FDR is much larger than the imposed adjusted p-value threshold, the fraction of false discoveries is not controlled at the claimed level. \n", sep = ''), resultfile)
    writeLines(c(paste("```{r fdr, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'fdr')",
                 "```", "---"), resultfile)
  }
  
  if ('fdrvsexpr' %in% setup.parameters$comparisons) {
    writeLines("<a name='fdrvsexpr'></a>", resultfile)
    writeLines("## FDR vs average expression level [(Contents)](#contents)", resultfile)
    writeLines(paste("The figures below show the observed false discovery rate as a function of the the average expression level of the genes ('A', shown on the x-axis). The average expression levels in a given data set are binned into 10 quantiles (each containing 1/10 of the values) and the false discovery rate at an imposed adjusted p-value cutoff of ", setup.parameters$fdr.threshold, "is estimated for each bin. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n"), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"),
                   resultfile)
      }
      writeLines(c(paste("```{r fdrvsexpr-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ",
                         6*nbr.cols2, ", fig.height = ", 6*nbr.rows2, ", fig.align = 'left'}", sep = ""),
                   paste("plotFDRVsExpr(setup.parameters, sel.nbrsamples = ",
                         setup.parameters$incl.nbr.samples[i], ", sel.repl =  c(", paste(intersect(setup.parameters$incl.replicates, setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)]), collapse = ", "), "))", sep = ""),
                   "```", "---"), resultfile)
    }
  }
  
  if ('tpr' %in% setup.parameters$comparisons) {
    writeLines("<a name='tpr'></a>", resultfile)
    writeLines("## TPR [(Contents)](#contents)", resultfile)
    writeLines(paste("The true positive rate (TPR) indicates the fraction of truly non-differentially expressed genes that are indeed considered significant by a method at a given significance threshold. A good method gives a high true positive rate, while at the same time keeping the false discovery rate under control. The figures below show the observed rate of true positives at an adjusted p-value threshold of ", setup.parameters$tpr.threshold,". Only methods returning corrected p-values or FDR estimates are included. Each boxplot summarizes the values obtained across all data set replicates included in the comparison.\n", sep = ''), resultfile)
    writeLines(c(paste("```{r tpr, dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 14, fig.height =", length(setup.parameters$incl.nbr.samples) + 0.65*kk1, ", fig.align = 'left'}"),
                 "plotResultTable(setup.parameters, res.table, 'tpr')",
                 "```", "---"), resultfile)
  }
  
  if ('overlap' %in% setup.parameters$comparisons) {
    writeLines("<a name='overlap'></a>", resultfile)
    writeLines("## Overlap, single replicate [(Contents)](#contents)", resultfile)
    writeLines(paste("The table below shows, for each pair of differential expression methods, the number of genes that are considered statistically significant by both of them at an adjusted p-value threshold of ",setup.parameters$overlap.threshold,". Only methods returning corrected p-values or FDR estimates are included in the comparison. Note that the size of the overlap between two sets naturally depends on the number of genes in each of the sets (indicated along the diagonal of the table).\n", sep = ''), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        idx <- which(setup.parameters$file.info$nbr.samples == setup.parameters$incl.nbr.samples[i])
      } else {
        idx <- which(is.na(setup.parameters$file.info$nbr.samples))
      }
      if (length(unique(setup.parameters$file.info$de.methods[idx])) > 1) {
        if (any(!is.na(setup.parameters$incl.nbr.samples))) {
          writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"), resultfile)
        }
        writeLines(c(paste("```{r overlap-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = 10, fig.height = 10, fig.align = 'left'}", sep = ""),
##                     "options(width = 200)",
                     paste("computeOverlap(setup.parameters, sel.nbrsamples = ", setup.parameters$incl.nbr.samples[i], ", sel.repl = ", mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ", plot.type = 'Overlap')"),
                     "```", "---"), resultfile)
      }
    }
  }
  
  if ('sorensen' %in% setup.parameters$comparisons) {
    writeLines("<a name='sorensen'></a>", resultfile)
    writeLines("## Sorensen index, single replicate [(Contents)](#contents)", resultfile)
    writeLines(paste("The table below shows, for each pair of differential expression methods, the Sorensen index, which is a way of quantifying the overlap between the collections of differentially expressed genes found by the two methods at an adjusted p-value threshold of ", setup.parameters$overlap.threshold,".Only methods returning corrected p-values or FDR estimates are included. The Sorensen index is defined as the ratio between the number of genes shared by the two sets and the average number of genes in the two sets. Hence, it always attains values between 0 and 1. A larger Sorensen index implies better overlap between the two sets, and hence that the two compared methods give similar differential expression results. The values of the Sorensen index for all pairs of compared methods are also visualized in a 'heatmap', where the color corresponds to Sorensen index.\n", sep = ''), resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        idx <- which(setup.parameters$file.info$nbr.samples == setup.parameters$incl.nbr.samples[i])
      } else {
        idx <- which(is.na(setup.parameters$file.info$nbr.samples))
      }
      if (length(unique(setup.parameters$file.info$de.methods[idx])) > 1) {
        if (any(!is.na(setup.parameters$incl.nbr.samples))) {
          writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"), resultfile)
        }
        writeLines(c(paste("```{r sorensen-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ", length(setup.parameters$incl.de.methods) + 3, ", fig.height = ", length(setup.parameters$incl.de.methods) + 3, ", fig.align = 'left'}", sep = ""),
                     paste("computeOverlap(setup.parameters, sel.nbrsamples = ", setup.parameters$incl.nbr.samples[i], ", sel.repl = ", mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, y = setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ", plot.type = 'Sorensen')"),
                     "```", "---"), resultfile)
      }
    }
  }
  
  if ('correlation' %in% setup.parameters$comparisons) {
    writeLines("<a name='correlation'></a>", resultfile)
    writeLines("## Spearman correlation between scores [(Contents)](#contents)", resultfile)
    writeLines("The table below shows, for each pair of compared differential expression methods, the Spearman correlation between the scores that they assign to the genes. The value of the correlation is always between -1 and 1, and a high positive value of the Spearman correlation indicates that the compared methods rank the genes in a similar fashion. The results are also shown in a 'heatmap', where the color indicates the Spearman correlation. Finally, the methods are clustered using hierarchical clustering, with a dissimilarity measure defined as 1 - Spearman correlation. This visualizes the relationships among the compared differential expression methods, and groups together methods that rank the genes similarly.\n", resultfile)
    for (i in seq_len(length(setup.parameters$incl.nbr.samples))) {
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        idx <- which(setup.parameters$file.info$nbr.samples == setup.parameters$incl.nbr.samples[i])
      } else {
        idx <- which(is.na(setup.parameters$file.info$nbr.samples))
      }
      if (any(!is.na(setup.parameters$incl.nbr.samples))) {
        writeLines(paste("####", setup.parameters$incl.nbr.samples[i], " samples/condition [(Contents)](#contents)"), resultfile)
      }
      writeLines(c(paste("```{r correlation-", i, ", dev = c('png', 'pdf'), eval = TRUE, include = TRUE, fig.width = ", length(setup.parameters$incl.de.methods) + 3, ", fig.height = ", length(setup.parameters$incl.de.methods) + 3, ", fig.align = 'left'}", sep = ""),
                   paste("computeCorrelation(setup.parameters, sel.nbrsamples = ", setup.parameters$incl.nbr.samples[i], ", sel.repl = ", mfv(setup.parameters$file.info$repl[vapply(setup.parameters$file.info$nbr.samples, identical, setup.parameters$incl.nbr.samples[i], FUN.VALUE = FALSE)])[1], ")"),
                   "```", "---"), resultfile)
    }
  }
  
  close(resultfile)
}

checkClass <- function(object, objname, trueclass) {
  if (!(is(object, trueclass))) {
    stop(paste("The object", objname, "should be of class", trueclass))
  }
}

#' A GUI to the main function for running the performance comparison between differential expression methods. 
#' 
#' This function provides a GUI to the main function for performing comparisons among differential expression methods and generating a report in HTML format (\code{\link{runComparison}}). It is assumed that all differential expression results have been generated in advance (using e.g. the function \code{\link{runDiffExp}}) and that the result \code{compData} object for each data set and each differential expression method is saved separately in files with the extension \code{.rds}. The function opens a graphical user interface where the user can set parameter values and choose the files to be used as the basis of the comparison. It is, however, possible to circumvent the GUI and call the comparison function \code{\link{runComparison}} directly. 
#' 
#' This function requires that the \code{rpanel} package is installed. If this package can not be installed, please use the \code{\link{runComparison}} function directly.
#' 
#' @param input.directories A list of directories containing the result files (\code{*.rds}). All results in the provided directories will be available for inclusion in the comparison, and the selection is performed through a graphical user interface. All result objects saved in the files should be of the \code{compData} class, although list objects created by earlier versions of \code{compcodeR} are supported.
#' @param output.directory The directory where the results should be written. The subdirectory structure will be created automatically. If the directory already exists, it will be overwritten.
#' @param recursive A logical parameter indicating whether or not the search should be extended recursively to subfolders of the \code{input.directories}. 
#' @param out.width The width of the figures in the final report. Will be passed on to \code{knitr} when the HTML is generated. Can be for example "800px" (see \code{knitr} documentation for more information)
#' @param upper.limits,lower.limits Lists that can be used to manually set upper and lower limits for boxplots of fdr, tpr, auc, mcc, fracsign, nbrtpfp, nbrsign and typeIerror.
#' @return
#' The function will create a comparison report, named \strong{compcodeR_report<timestamp>.html}, in the \code{output.directory}. It will also create subfolders named \code{compcodeR_code} and \code{compcodeR_figure}, where the code used to perform the differential expression analysis and the figures contained in the report, respectively, will be saved. Note that if these directories already exist they will be overwritten.
#' @export
#' @author Charlotte Soneson
#' @examples
#' if (interactive()) {
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 12500, 
#'                                     samples.per.cond = 5, n.diffexp = 1250, 
#'                                     output.file = "mydata.rds")
#' runDiffExp(data.file = "mydata.rds", result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", output.directory = ".", 
#'            norm.method = "TMM")
#' runDiffExp(data.file = "mydata.rds", result.extent = "ttest", 
#'            Rmdfunction = "ttest.createRmd", output.directory = ".", 
#'            norm.method = "TMM")
#' runComparisonGUI(input.directories = ".", output.directory = ".", recursive = FALSE)
#' }
#nocov start
runComparisonGUI <- function(input.directories, output.directory, recursive, 
                             out.width = NULL, upper.limits = NULL, lower.limits = NULL) {
  if (!requireNamespace("rpanel", quietly = TRUE)) stop("To use the GUI you must install the rpanel package. If it is not possible to install the package on your system, please use the runComparison() function instead.")
  checkClass(input.directories, "input.directories", "character")
  checkClass(output.directory, "output.directory", "character")
  checkClass(recursive, "recursive", "logical")
  
  ## Get the absolute path of the input directories
  if (is.null(input.directories)) {
    stop("You have to provide at least one input directory!")
  } else {
    input.directories <- normalizePath(input.directories, winslash = "/")
  }
  
  if (is.null(output.directory)) {
    stop("You have to provide an output directory!")
  } else if (length(output.directory) != 1) {
    stop("Please provide precisely one output directory.")
  } else {
    output.directory <- normalizePath(output.directory, winslash = "/")
  }
  
  ## List all rds files in the input directories
  input.files.temp <- NULL
  for (input.directory in input.directories) {
    input.files.temp <- union(input.files.temp, file.path(input.directory,
                                                          list.files(input.directory,
                                                                     recursive = recursive,
                                                                     pattern = "*.rds$")))
  }
  input.files.temp <- normalizePath(input.files.temp, winslash = "/")
  if (length(input.files.temp) == 0) {
    stop("No .rds files in the input directory")
  }

  ## Go through the files in the input.directories and see which
  ## datasets are available.
  suppressWarnings({
    datasets <- NULL
    input.files <- NULL
    for (input.file in input.files.temp) {
      temp.input <- readRDS(input.file)
      ## Check if compatible list, then convert to compData object
      if (is.list(temp.input)) {
        temp.input <- convertListTocompData(temp.input)
      }
      ## Keep only the result files (no data files or other things).
      if (is.logical(check_compData_results(temp.input))) {
        datasets <- c(datasets, info.parameters(temp.input)$dataset) 
        input.files <- c(input.files, input.file)
      }
    }
    avail.datasets <- c('', unique(datasets))
  })
  
  if (length(avail.datasets) == 1) {
    stop("No result .rds files in the input directory")
  }
  
  ## Create the GUI for selecting the data set
  dataset.panel <- rpanel::rp.control(title = "Select data set", 
                                      panelname = "dataset.panel", 
                                      incl.dataset = avail.datasets[1],
                                      fdr.threshold = 0.05,
                                      tpr.threshold = 0.05,
                                      mcc.threshold = 0.05,
                                      typeI.threshold = 0.05,
                                      fdc.maxvar = 1500,
                                      overlap.threshold = 0.05,
                                      fracsign.threshold = 0.05,
                                      nbrtpfp.threshold = 0.05,
                                      ma.threshold = 0.05, 
                                      signal.measure = "mean",
                                      input.files = input.files,
                                      datasets = datasets, 
                                      avail.datasets = avail.datasets, 
                                      output.directory = output.directory, 
                                      out.width = out.width, 
                                      upper.limits = upper.limits, 
                                      lower.limits = lower.limits)
  
  ## Create the radiobutton for selecting a data set
  combo.dataset <- rpanel::rp.combo(dataset.panel,
                                    variable = "incl.dataset",
                                    action = I,
                                    vals = avail.datasets,
                                    initval = avail.datasets[1],
                                    prompt = "Data set",
                                    pos = list("row" = 1, "column" = 1))
  
  go.button <- rpanel::rp.button(dataset.panel, action = function(panel) {
    ## Go through all input files corresponding to the chosen data set, and record the 
    ## available nbrsamples/repl/de methods
    panel$input.files <- panel$input.files[panel$datasets == panel$incl.dataset]
    panel$datasets <- panel$datasets[panel$datasets == panel$incl.dataset]
    suppressWarnings({
      nbr.samples <- NULL
      repl <- NULL
      de.methods <- NULL
      for (input.file in panel$input.files) {
        temp.input <- readRDS(input.file)
        if (is.list(temp.input)) temp.input <- convertListTocompData(temp.input)
        nbr.samples <- c(nbr.samples, ifelse(is.null(info.parameters(temp.input)$samples.per.cond), 
                                             NA_real_, info.parameters(temp.input)$samples.per.cond))
        repl <- c(repl, ifelse(is.null(info.parameters(temp.input)$repl.id), 
                               NA_real_, info.parameters(temp.input)$repl.id))
        de.methods <- c(de.methods, ifelse(is.null(method.names(temp.input)$full.name), 
                                           NA_character_, method.names(temp.input)$full.name))
      }
      avail.nbr.samples <- sort(as.numeric(unique(nbr.samples)))  # numeric(0) if all are NA
      avail.repl <- sort(as.numeric(unique(repl)))  # numeric(0) if all are NA
      avail.de.methods <- setdiff(unique(de.methods), c(NA_character_))
    })
    panel$avail.nbr.samples <- avail.nbr.samples
    panel$avail.repl <- avail.repl
    panel$avail.de.methods <- avail.de.methods
    panel$nbr.samples <- nbr.samples
    panel$repl <- repl
    panel$de.methods <- de.methods
    createSelectionPanel(panel)
    return(panel)
  }, quitbutton = TRUE, title = "Continue", pos = list("row" = 2, "column" = 2))
}
#nocov end

# Create the user interface for selecting parameters for the method comparison
# 
# Launch the user interface for selecting the parameters that govern the comparison of differential expression methods. This function will be called when pressing the "Continue" button in the interface for selecting the data set. This function will not be called directly by the user. 
# 
# @param panel An rpanel
# @author Charlotte Soneson
#nocov start
createSelectionPanel = function(panel) {
  ###########################################
  ## Create the GUI for setting parameters ##
  ###########################################
  
  ## Define the set of possible comparison methods
  comparisons <- c("AUC", "ROC, one replicate", "ROC, all replicates",
                   "Type I error", "FDR", 
                   "FDR vs average expression level",
                   "TPR", "MCC", "Number significant", 
                   "Number of TP, FP, TN, FN",
                   "False discovery curves, one replicate",
                   "False discovery curves, all replicates",
                   "Fraction significant", "Overlap, one replicate",
                   "MA plot", "Gene score vs average expression level",
                   "Gene score vs signal for condition-specific genes",
                   "Score distribution vs number of outliers", 
                   "Spearman correlation between scores",
                   "Sorensen index, one replicate")
  
  ## Create the main panel window
  if (length(panel$avail.nbr.samples) != 0) {
    incl.nbr.samples <- rep(TRUE, length(panel$avail.nbr.samples))
  } else {
    incl.nbr.samples <- NULL
  }
  if (length(panel$avail.repl) != 0) {
    incl.repl <- rep(TRUE, length(panel$avail.repl))
  } else {
    incl.repl <- NULL
  }
  main.panel <- rpanel::rp.control(title = "Set parameters",
                                   panelname = "main.panel",
                                   incl.nbr.samples = incl.nbr.samples,
                                   incl.replicates = incl.repl,
                                   incl.de.methods = rep(TRUE, length(panel$avail.de.methods)),
                                   incl.dataset = panel$incl.dataset,
                                   fdr.threshold = panel$fdr.threshold,
                                   tpr.threshold = panel$tpr.threshold,
                                   mcc.threshold = panel$mcc.threshold,
                                   typeI.threshold = panel$typeI.threshold,
                                   fdc.maxvar = panel$fdc.maxvar,
                                   ma.threshold = panel$ma.threshold,
                                   signal.measure = panel$signal.measure,
                                   overlap.threshold = panel$overlap.threshold,
                                   fracsign.threshold = panel$fracsign.threshold,
                                   nbrtpfp.threshold = panel$nbrtpfp.threshold,
                                   comparisons = rep(TRUE, length(comparisons)),
                                   input.files = panel$input.files,
                                   nbr.samples = panel$nbr.samples,
                                   repl = panel$repl,
                                   de.methods = panel$de.methods,
                                   datasets = panel$datasets, 
                                   output.directory = panel$output.directory, 
                                   out.width = panel$out.width, 
                                   upper.limits = panel$upper.limits, 
                                   lower.limits = panel$lower.limits)
  
  ## Create the checkboxes for selecting the number of samples
  if (length(panel$avail.nbr.samples) != 0) {
    checkbox.incl.nbr.samples <- 
      rpanel::rp.checkbox(main.panel,
                          variable = "incl.nbr.samples",
                          action = I,
                          labels = panel$avail.nbr.samples,
                          initval = labels %in% main.panel$incl.nbr.samples,
                          title = "Number of samples",
                          pos = list("row" = 3, "column" = 2, 
                                     "rowspan" = length(panel$avail.nbr.samples)))
  }
  
  ## Create the checkboxes for selecting the replicates to include
  if (length(panel$avail.repl) != 0) {
    checkbox.repl <- 
      rpanel::rp.checkbox(main.panel,
                          variable = "incl.replicates",
                          action = I,
                          labels = panel$avail.repl,
                          initval = labels %in% main.panel$incl.replicates,
                          title = "Replicates",
                          pos = list("row" = 2, "column" = 2))
  }
  
  ## Create a checkbox for selecting the DE methods to include
  checkbox.de.methods <- 
    rpanel::rp.checkbox(main.panel,
                        variable = "incl.de.methods",
                        action = I,
                        labels = panel$avail.de.methods,
                        initval = labels %in% main.panel$incl.de.methods,
                        title = "DE methods",
                        pos = list("row" = 2, "column" = 1))
  
  ## Create a checkbox for selecting the comparison methods to include
  checkbox.comparisons <- 
    rpanel::rp.checkbox(main.panel,
                        variable = "comparisons",
                        action = I,
                        labels = comparisons,
                        initval = labels %in% main.panel$comparisons,
                        title = "Comparisons",
                        pos = list("row" = 2, "column" = 3))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## used for the FDR analysis
  text.fdr.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "fdr.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for FDR",
                         initval = main.panel$fdr.threshold,
                         title = "",
                         pos = list("row" = 3, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## used for the TPR analysis
  text.tpr.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "tpr.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for TPR",
                         initval = main.panel$tpr.threshold,
                         title = "",
                         pos = list("row" = 4, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## used for the MCC analysis
  text.mcc.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "mcc.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for MCC",
                         initval = main.panel$mcc.threshold,
                         title = "",
                         pos = list("row" = 4, "column" = 1))
  
  ## Create a text field for entering the nominal p-value threshold
  ## used for the type I error analysis
  text.typeI.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "typeI.threshold",
                         action = I,
                         labels = "Nominal p-value threshold for type I error",
                         initval = main.panel$typeI.threshold,
                         title = "",
                         pos = list("row" = 5, "column" = 1))
  
  ## Create a text field for entering the maximal number of variables
  ## to use in the false discovery curves
  text.fdc.maxvar <- 
    rpanel::rp.textentry(main.panel,
                         variable = "fdc.maxvar",
                         action = I,
                         labels = "Maximal number of top variables in false discovery curves",
                         initval = main.panel$fdc.maxvar,
                         title = "",
                         pos = list("row" = 6, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## for overlap analysis
  text.overlap.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "overlap.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for overlap analysis",
                         initval = main.panel$overlap.threshold,
                         title = "",
                         pos = list("row" = 7, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## for analysis of the fraction of significant variables
  text.fracsign.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "fracsign.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for analysis of fraction significant variables",
                         initval = main.panel$fracsign.threshold,
                         title = "",
                         pos = list("row" = 8, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## for analysis of the number of TP, FP, TN, FN
  text.nbrtpfp.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "nbrtpfp.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for analysis of the number of TP, FP, TN, FN variables",
                         initval = main.panel$nbrtpfp.threshold,
                         title = "",
                         pos = list("row" = 9, "column" = 1))
  
  ## Create a text field for entering the adjusted p-value threshold
  ## for coloring MA plots
  text.ma.threshold <- 
    rpanel::rp.textentry(main.panel,
                         variable = "ma.threshold",
                         action = I,
                         labels = "Adjusted p-value threshold for coloring of genes in MA plots",
                         initval = main.panel$ma.threshold,
                         title = "",
                         pos = list("row" = 10, "column" = 1))
  
  ## Create a radiobutton to select the signal measure 
  combo.signal.measure <- rpanel::rp.combo(main.panel,
                                           variable = "signal.measure",
                                           action = I,
                                           vals = c("mean", "snr"),
                                           initval = main.panel$signal.measure,
                                           prompt = "Signal measure for condition-specific genes",
                                           pos = list("row" = 11, "column" = 1))
  
  ## Create the "OK" button that runs the comparison
  ok.button <- 
    rpanel::rp.button(main.panel, action = performComparison,
                      title = "Run comparison",
                      quitbutton = TRUE,
                      pos = list("row" = 3, "column" = 3, "rowspan" = 1))
  
  ## Create the "Close" button
  close.button <- 
    rpanel::rp.button(main.panel, action = function(panel){return(panel)}, 
                      title = "Cancel", quitbutton = TRUE, 
                      pos = list("row" = 4, "column" = 3, "rowspan" = 1))
  
  #rp.block(main.panel)
}
#nocov end

checkRange <- function(value, parname, minvalue, maxvalue) {
  if (is.na(as.numeric(value))) {
    stop(paste("Illegal value for parameter", parname))
  } else {
    newvalue <- as.numeric(value)
    if (newvalue < minvalue) newvalue <- minvalue
    if (newvalue > maxvalue) newvalue <- maxvalue
  }
  return(newvalue)
}

# Run the comparison of differential expression methods
# 
# This function will be launched when pressing the "Run comparison" button in the parameter selection user interface. It will not be called directly by the user. 
# 
# @param panel An rpanel
# @author Charlotte Soneson
#nocov start
performComparison <- function(panel) {
  message("Be patient, your analysis is running...")
  
  ## Extract all parameters from the panel object
  suppressWarnings({
    if (length(panel$incl.nbr.samples) == 1) {
      names(panel$incl.nbr.samples) <- 
        unique(panel$nbr.samples[!is.na(panel$nbr.samples)])
    }
    if (length(panel$incl.replicates) == 1) {
      names(panel$incl.replicates) <- unique(panel$repl[!is.na(panel$repl)])
    }
    if (length(panel$incl.de.methods) == 1) {
      names(panel$incl.de.methods) <- unique(panel$de.methods)
    }
    parameters <- list(incl.nbr.samples = as.numeric(names(which(panel$incl.nbr.samples == TRUE))),
                       incl.replicates = as.numeric(names(which(panel$incl.replicates == TRUE))),
                       incl.dataset = panel$incl.dataset,
                       incl.de.methods = names(which(panel$incl.de.methods == TRUE)),
                       fdr.threshold = as.numeric(panel$fdr.threshold),
                       tpr.threshold = as.numeric(panel$tpr.threshold),
                       mcc.threshold = as.numeric(panel$mcc.threshold),
                       typeI.threshold = as.numeric(panel$typeI.threshold),
                       ma.threshold = as.numeric(panel$ma.threshold), 
                       fdc.maxvar = as.numeric(panel$fdc.maxvar),
                       signal.measure = as.character(panel$signal.measure),
                       overlap.threshold = as.numeric(panel$overlap.threshold),
                       fracsign.threshold = as.numeric(panel$fracsign.threshold),
                       nbrtpfp.threshold = as.numeric(panel$nbrtpfp.threshold),
                       upper.limits = panel$upper.limits, 
                       lower.limits = panel$lower.limits,
                       comparisons = names(which(panel$comparisons == TRUE)))
    nbr.samples <- as.numeric(panel$nbr.samples)
    repl <- as.numeric(panel$repl)
    de.methods <- panel$de.methods
    datasets <- panel$datasets
    input.files <- panel$input.files
    output.directory <- panel$output.directory
    out.width <- panel$out.width
  })
  
  ## Check ranges etc for thresholds
  parameters$fdr.threshold <- checkRange(parameters$fdr.threshold, "fdr.treshold", 0, 1)
  parameters$tpr.threshold <- checkRange(parameters$tpr.threshold, "tpr.threshold", 0, 1)
  parameters$mcc.threshold <- checkRange(parameters$mcc.threshold, "mcc.threshold", 0, 1)
  parameters$typeI.threshold <- checkRange(parameters$typeI.threshold, "typeI.threshold", 0, 1)
  parameters$ma.threshold <- checkRange(parameters$ma.threshold, "ma.threshold", 0, 1)
  parameters$overlap.threshold <- checkRange(parameters$overlap.threshold, "overlap.threshold", 0, 1)
  parameters$fracsign.threshold <- checkRange(parameters$fracsign.threshold, "fracsign.threshold", 0, 1)
  parameters$nbrtpfp.threshold <- checkRange(parameters$nbrtpfp.threshold, "nbrtpfp.threshold", 0, 1)
  
  ## Transform the names of the comparisons to make
  parameters$comparisons <- shorten.method.names(parameters$comparisons)
  
  ## Run the comparison 
  file.table <- data.frame(input.files = input.files, 
                           datasets = datasets,
                           nbr.samples = nbr.samples, 
                           repl = repl, 
                           de.methods = de.methods, 
                           stringsAsFactors = FALSE)
  runComparison(file.table, parameters, 
                output.directory, check.table = FALSE, 
                out.width = out.width)
  
  return(panel)
}
#nocov end

#' Run the performance comparison between differential expression methods. 
#' 
#' The main function for performing comparisons among differential expression methods and generating a report in HTML format. It is assumed that all differential expression results have been generated in advance (using e.g. the function \code{\link{runDiffExp}}) and that the result \code{compData} object for each data set and each differential expression method is saved separately in files with the extension \code{.rds}. Note that the function can also be called via the \code{\link{runComparisonGUI}} function, which lets the user set parameters and select input files using a graphical user interface. 
#' 
#' The input to \code{\link{runComparison}} is a data frame with at least a column named \code{input.files}, containing paths to \code{.rds} files containing result objects (of the class \code{compData}), such as those generated by \code{\link{runDiffExp}}. Other columns that can be included in the data frame are \code{datasets}, \code{nbr.samples}, \code{repl} and \code{de.methods}. They have to match the information contained in the corresponding result objects. If these columns are not present, they will be added to the data frame automatically.
#' 
#' @param file.table A data frame with at least a column \code{input.files}, potentially also columns named \code{datasets}, \code{nbr.samples}, \code{repl} and \code{de.methods}.
#' @param parameters A list containing parameters for the comparison study. The following entries are supported, and used by different comparison methods:
#' \itemize{
#' \item \code{incl.nbr.samples} An array with sample sizes (number of samples per condition) to consider in the comparison. If set to \code{NULL}, all sample sizes will be included.
#' \item \code{incl.dataset} A dataset name (corresponding to the \code{dataset} slot of the results or data objects), indicating the dataset that will be used for the comparison. Only one dataset can be chosen. 
#' \item \code{incl.replicates} An array with replicate numbers to consider in the comparison. If set to \code{NULL}, all replicates will be included.
#' \item \code{incl.de.methods} An array with differential expression methods to be compared. If set to \code{NULL}, all differential expression methods will be included.
#' \item \code{fdr.threshold} The adjusted p-value threshold for FDR calculations. Default 0.05.
#' \item \code{tpr.threshold} The adjusted p-value threshold for TPR calculations. Default 0.05.
#' \item \code{mcc.threshold} The adjusted p-value threshold for MCC calculations. Default 0.05.
#' \item \code{typeI.threshold} The nominal p-value threshold for type I error calculations. Default 0.05.
#' \item \code{fdc.maxvar} The maximal number of variables to include in false discovery curve plots. Default 1500.
#' \item \code{overlap.threshold} The adjusted p-value for overlap analysis. Default 0.05.
#' \item \code{fracsign.threshold} The adjusted p-value for calculation of the fraction/number of genes called significant. Default 0.05.
#' \item \code{nbrtpfp.threshold} The adjusted p-value for calculation of the number of TP, FP, TN, FN genes. Default 0.05.
#' \item \code{ma.threshold} The adjusted p-value threshold for coloring genes in MA plots. Default 0.05.
#' \item \code{signal.measure} Either \code{'mean'} or \code{'snr'}, determining how to define the signal strength for a gene which is expressed in only one condition.
#' \item \code{upper.limits,lower.limits} Lists that can be used to manually set the upper and lower plot limits for boxplots of fdr, tpr, auc, mcc, fracsign, nbrtpfp and typeIerror.
#' \item \code{comparisons} Array containing the comparison methods to be applied. The entries must be chosen among the following abbreviations:
#' \itemize{
#' \item \code{"auc"} - Compute the area under the ROC curve
#' \item \code{"mcc"} - Compute Matthew's correlation coefficient
#' \item \code{"tpr"} - Compute the true positive rate at a given adjusted p-value threshold (\code{tpr.threshold}) 
#' \item \code{"fdr"} - Compute the false discovery rate at a given adjusted p-value threshold (\code{fdr.threshold}) 
#' \item \code{"fdrvsexpr"} - Compute the false discovery rate as a function of the expression level.
#' \item \code{"typeIerror"} - Compute the type I error rate at a given nominal p-value threshold (\code{typeI.threshold}) 
#' \item \code{"fracsign"} - Compute the fraction of genes called significant at a given adjusted p-value threshold (\code{fracsign.threshold}).
#' \item \code{"nbrsign"} - Compute the number of genes called significant at a given adjusted p-value threshold (\code{fracsign.threshold}).
#' \item \code{"nbrtpfp"} - Compute the number of true positives, false positives, true negatives and false negatives at a given adjusted p-value threshold (\code{nbrtpfp.threshold}).
#' \item \code{"maplot"} - Construct MA plots, depicting the average expression level and the log fold change for the genes and indicating the genes called differential expressed at a given adjusted p-value threshold (\code{ma.threshold}).
#' \item \code{"fdcurvesall"} - Construct false discovery curves for each of the included replicates.
#' \item \code{"fdcurvesone"} - Construct false discovery curves for a single replicate only
#' \item \code{"rocall"} - Construct ROC curves for each of the included replicates 
#' \item \code{"rocone"} - Construct ROC curves for a single replicate only
#' \item \code{"overlap"} - Compute the overlap between collections of genes called differentially expressed by the different methods at a given adjusted p-value threshold (\code{overlap.threshold})
#' \item \code{"sorensen"} - Compute the Sorensen index, quantifying the overlap between collections of genes called differentially expressed by the different methods, at a given adjusted p-value threshold (\code{overlap.threshold})
#' \item \code{"correlation"} - Compute the Spearman correlation between gene scores assigned by different methods 
#' \item \code{"scorevsoutlier"} - Visualize the distribution of the gene scores as a function of the number of outlier counts introduced for the genes
#' \item \code{"scorevsexpr"} - Visualize the gene scores as a function of the average expression level of the genes
#' \item \code{"scorevssignal"} - Visualize the gene score as a function of the 'signal strength' (see the \code{signal.measure} parameter above) for genes that are expressed in only one condition
#' }
#' }
#' @param output.directory The directory where the results should be written. The subdirectory structure will be created automatically. If the directory already exists, it will be overwritten.
#' @param check.table Logical, should the input table be checked for consistency. Default \code{TRUE}.
#' @param out.width The width of the figures in the final report. Will be passed on to \code{knitr} when the HTML is generated.
#' @param saveResultTable Logical, should the intermediate result table be saved for future use ? Default to \code{FALSE}.
#' @param knitResults Logical, should the Rmd be generated and knited ? Default to \code{TRUE}.
#' @return
#' The function will create a comparison report, named \strong{compcodeR_report<timestamp>.html}, in the \code{output.directory}. It will also create subfolders named \code{compcodeR_code} and \code{compcodeR_figure}, where the code used to perform the differential expression analysis and the figures contained in the report, respectively, will be stored. Note that if these directories already exists, they will be overwritten.
#' @export
#' @author Charlotte Soneson
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", output.directory = tmpdir, 
#'            norm.method = "TMM")
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.exact", 
#'            Rmdfunction = "edgeR.exact.createRmd", output.directory = tmpdir, 
#'            norm.method = "TMM", 
#'            trend.method = "movingave", disp.type = "tagwise")
#' file.table <- data.frame(input.files = file.path(tmpdir, 
#'                          c("mydata_voom.limma.rds", "mydata_edgeR.exact.rds")), 
#'                          stringsAsFactors = FALSE)
#' parameters <- list(incl.nbr.samples = 5, incl.replicates = 1, incl.dataset = "mydata", 
#'                    incl.de.methods = NULL, 
#'                    fdr.threshold = 0.05, tpr.threshold = 0.05, typeI.threshold = 0.05,
#'                    ma.threshold = 0.05, fdc.maxvar = 1500, overlap.threshold = 0.05,
#'                    fracsign.threshold = 0.05, mcc.threshold = 0.05, 
#'                    nbrtpfp.threshold = 0.05, 
#'                    comparisons = c("auc", "fdr", "tpr", "ma", "correlation"))
#' if (interactive()) {
#'   runComparison(file.table = file.table, parameters = parameters, output.directory = tmpdir)
#' }
#' @importFrom grDevices heat.colors
#' @importFrom graphics axis legend lines par title
#' @importFrom stats as.dist cor hclust loess median na.omit predict rexp rnbinom
#'   rpois runif sd
#' @importFrom utils packageVersion
runComparison <- function(file.table, 
                          parameters, 
                          output.directory, 
                          check.table = TRUE, 
                          out.width = NULL,
                          saveResultTable = FALSE,
                          knitResults = TRUE) {
  
  checkClass(file.table, "file.table", "data.frame")
  checkClass(output.directory, "output.directory", "character")
  
  if (is.null(parameters)) {
    parameters <- list()
  } 
  checkClass(parameters, "parameters", "list")
  
  ## Set default values if not given
  if (is.null(parameters$fdr.threshold)) parameters$fdr.threshold <- 0.05
  if (is.null(parameters$tpr.threshold)) parameters$tpr.threshold <- 0.05
  if (is.null(parameters$mcc.threshold)) parameters$mcc.threshold <- 0.05
  if (is.null(parameters$typeI.threshold)) parameters$typeI.threshold <- 0.05
  if (is.null(parameters$ma.threshold)) parameters$ma.threshold <- 0.05
  if (is.null(parameters$fdc.maxvar)) parameters$fdc.maxvar <- 1500
  if (is.null(parameters$signal.measure)) parameters$signal.measure <- "mean"
  if (is.null(parameters$overlap.threshold)) parameters$overlap.threshold <- 0.05
  if (is.null(parameters$fracsign.threshold)) parameters$fracsign.threshold <- 0.05
  if (is.null(parameters$nbrtpfp.threshold)) parameters$nbrtpfp.threshold <- 0.05
  if (is.null(parameters$comparisons)) parameters$comparisons <- 
    c("auc", "fdr", "tpr", "mcc", 
      "maplot", "correlation", "nbrtpfp", 
      "typeIerror", "fracsign", "nbrsign", 
      "fdcurvesall", "fdcurvesone", 
      "rocall", "rocone", "overlap", 
      "scorevsexpr", "sorensen", 
      "correlation", "scorevsoutlier",
      "fdrvsexpr", "scorevssignal")
  
  ## Some options, useful for debugging
  save.rmdfile <- FALSE  ## set to TRUE to save the .Rmd file
  
  file.table <- data.frame(file.table, stringsAsFactors = FALSE)
  
  file.table$input.files <- normalizePath(as.character(file.table$input.files),
                                          winslash = "/")
  
  ## Check the input table/parameters
  ## First, make the table consistent (go through all input files, check the number of samples, 
  ## the replicate ID and the de.method. If it doesn't agree with table, change the table and 
  ## print a message.)
  ## This is not needed if we call runComparison from runComparisonGUI
  if (check.table) {
    file.table <- checkTableConsistency(file.table)
  }
  if (nrow(file.table) == 0) {
    stop("No methods left to compare after checking object validity!")
  }
  if (is.null(parameters$incl.dataset)) {
    parameters$incl.dataset <- unique(file.table$datasets)
  }
    
  ## Get the whole path to the output.directory
  output.directory <- normalizePath(output.directory, winslash = "/")
  
  ## If the output.directory does not exist, create it and some subfolders
  dir.create(output.directory, showWarnings = FALSE)
  dir.create(file.path(output.directory, "compcodeR_code"), showWarnings = FALSE)
  dir.create(file.path(output.directory, "compcodeR_figure"), showWarnings = FALSE)
  
  ## Set options for knitr
  opts_knit$set(progress = FALSE, verbose = FALSE)
  opts_chunk$set(warning = FALSE, error = FALSE, message = FALSE, comment = '',
                 fig.path = file.path(output.directory, "compcodeR_figure/"),
                 echo = FALSE, out.width = out.width)
  
  ## Go through the input files and keep only the ones that are compatible
  ## with the selected parameters. 
  ## First select the ones from the right dataset (this is a required parameter)
  idx.keep <- which(file.table$datasets %in% parameters$incl.dataset)
  file.table <- file.table[idx.keep, ]
  if (nrow(file.table) == 0) {
    stop("No methods left to compare after matching with datasets to include!")
  }
  
  ## If any of the files have information about nbr.samples, keep only the files
  ## with nbr.samples in accordance with the selected nbr.samples. Otherwise, 
  ## keep all files
  if (is.null(parameters$incl.nbr.samples) || all(is.na(parameters$incl.nbr.samples))) {
    ## If no selection of nbr.samples, keep all
    idx.keep <- seq_len(nrow(file.table))
  } else {
    if (any(!is.na(file.table$nbr.samples))) {
      idx.keep <- which(file.table$nbr.samples %in% parameters$incl.nbr.samples)
    } else {
      idx.keep <- seq_len(nrow(file.table))
    }
  }
  file.table <- file.table[idx.keep, ]
  if (nrow(file.table) == 0) {
    stop("No methods left to compare after matching with nbr.samples to include!")
  }
  
  ## If any of the files have information about replicate, keep only the files
  ## with replicate in accordance with the selected replicate. Otherwise, 
  ## keep all files
  if (is.null(parameters$incl.replicates) || all(is.na(parameters$incl.replicates))) {
    ## If no selection of replicates, keep all
    idx.keep <- seq_len(nrow(file.table))
  } else {
    if (any(!is.na(file.table$repl))) {
      idx.keep <- which(file.table$repl %in% parameters$incl.replicates)
    } else {
      idx.keep <- seq_len(nrow(file.table))
    }
  }
  file.table <- file.table[idx.keep, ]
  if (nrow(file.table) == 0) {
    stop("No methods left to compare after matching with replicates to include!")
  }
  
  if (is.null(parameters$incl.de.methods) || all(is.na(parameters$incl.de.methods))) {
    ## If no selection of de.methods, keep all
    idx.keep <- seq_len(nrow(file.table))
  } else {
    idx.keep <- which(file.table$de.methods %in% parameters$incl.de.methods)
  }
  file.table <- file.table[idx.keep, ]
  if (nrow(file.table) == 0) {
    stop("No methods left to compare after matching with DE methods to include!")
  }
  
  parameters$file.info <- file.table
  parameters$file.info[is.na(parameters$file.info)] <- NA_real_
  
  parameters$incl.nbr.samples <- unique(file.table$nbr.samples)
  parameters$incl.nbr.samples[is.na(parameters$incl.nbr.samples)] <- NA_real_
  parameters$incl.replicates <- unique(file.table$repl)
  parameters$incl.replicates[is.na(parameters$incl.replicates)] <- NA_real_
  parameters$incl.de.methods <- unique(file.table$de.methods)
  
  ## Check that the input files actually correspond to the same simulation settings
  ## (data sets with the same name can only differ in the number of samples and 
  ## the replicate number)
  checkCompatibility(parameters$file.info)
  
  parameters$specified.colors <- definePlotColors(parameters$incl.de.methods)
  
  setup.parameters <- parameters
  ##saveRDS(setup.parameters, "/Users/charlotte/Desktop/setup.parameters.ex.rds")

  ## Generate the Rmd file containing the results
  timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
  
  ## Save the setup.parameters
  setup.parameters.file <- file.path(output.directory, 
                                     paste("compcodeR_parameters_", 
                                           timestamp, 
                                           ".rds", sep = ""))
  saveRDS(setup.parameters, setup.parameters.file)
  
  if (saveResultTable) {
    doSaveResultsTable(setup.parameters.file,
                       output.file = file.path(output.directory, 
                                               paste("compcodeR_result_table_", 
                                                     timestamp, 
                                                     ".rds", sep = "")))
  }
  
  if (knitResults) {
    createResultsRmdFile(setup.parameters.file,
                         output.file = file.path(output.directory, 
                                                 paste("compcodeR_report_", 
                                                       timestamp, 
                                                       ".Rmd", sep = "")))
    
    ## Generate the HTML file
    out <- knitr::knit(input = file.path(output.directory, paste("compcodeR_report_",
                                                                 timestamp, 
                                                                 ".Rmd", sep = "")),
                       output = file.path(output.directory, paste("compcodeR_report_",
                                                                  timestamp, 
                                                                  ".md", sep = "")),
                       envir =  new.env())
    markdown::markdownToHTML(file = out,
                             output = file.path(output.directory, paste("compcodeR_report_",
                                                                        timestamp, 
                                                                        ".html", sep = "")),
                             encoding = "UTF-8",
                             title = "compcodeR results")
    
    ## Remove the .Rmd file
    if (!save.rmdfile) {
      file.remove(file.path(output.directory, paste("compcodeR_report_",
                                                    timestamp, 
                                                    ".Rmd", sep = "")))
      file.remove(file.path(output.directory, paste("compcodeR_report_",
                                                    timestamp, 
                                                    ".md", sep = "")))
    }
    
    ## Generate the HTML files containing the code for all compared instances
    generateCodeHTMLs(setup.parameters$file.info$input.files, 
                      file.path(output.directory, "compcodeR_code"))
  }
  message("Done!")
}

# Save the result table for further ploting
# 
# This function is generally not called by the user.
# 
# @author Charlotte Soneson
doSaveResultsTable <- function(setup.parameters.file, output.file) {
  ## Check that the output file ends with .rds
  if (!(substr(output.file, nchar(output.file) - 3, nchar(output.file)) == ".rds")) {
    output.file <- sub(strsplit(output.file, "\\.")[[1]][length(strsplit(output.file, "\\.")[[1]])], 
                       "rds", 
                       output.file)
  }
  
  ## Load the setup parameters
  setup.parameters <- readRDS(setup.parameters.file)
  setup.parameters$incl.nbr.samples <- as.numeric(setup.parameters$incl.nbr.samples)
  setup.parameters$incl.replicates <- as.numeric(setup.parameters$incl.replicates)
  if (all(!is.na(setup.parameters$incl.nbr.samples))) {
    setup.parameters$incl.nbr.samples <- sort(setup.parameters$incl.nbr.samples)
  }
  if (all(!is.na(setup.parameters$incl.replicates))) {
    setup.parameters$incl.replicates <- sort(setup.parameters$incl.replicates)
  }
  
  ## Create the result table
  if (any(c('auc', 'fdr', 'tpr', 'typeIerror', 'fracsign', 'nbrtpfp', 'nbrsign', 'mcc') %in% setup.parameters$comparisons)) {
    res.table <- createResultTable(setup.parameters)
    if (any(!is.na(setup.parameters$incl.nbr.samples))) {
      res.table <- padResultTable(res.table)
    }
  }
  
  ## Create the result table
  saveRDS(res.table, file = output.file)
}

# Check the compatibility of data sets used for method comparison 
# 
# Check that there are no contradictions or inconsistencies among the data sets used for comparing different methods for differential expression analysis. Please note that data sets that are given the same name (the 'dataset' entry) can only differ in the number of samples per condition and the replicate ID (and the unique data set ID). All other simulation parameters must be identical. This function is generally not called by the user.
# 
# @param file.info.table Table containing the files that are intended for use in the comparison. 
# @author Charlotte Soneson
checkCompatibility <- function(file.info.table){
  ## First check that all files have the same data set name
  if (length(unique(file.info.table$datasets[!is.na(file.info.table$datasets)])) != 1) {
    stop("All data sets must have the same name (the 'dataset' entry).")
  }
  
  ## Then check that all input files with the same nbr samples and repl.id are
  ## identical (that they have the same unique ID)
  all.combinations <- interaction(file.info.table$nbr.samples, 
                                  file.info.table$repl)
  for (tmp in levels(all.combinations)) {
    tmp.nbr.samples <- strsplit(tmp, "\\.")[[1]][1]
    tmp.repl <- strsplit(tmp, "\\.")[[1]][2]
    tmp.uID <- NULL
    tmp.lines <- which(all.combinations == tmp)
    if (length(tmp.lines) != 0) {
      for (i in seq_len(length(tmp.lines))) {
        w <- readRDS(as.character(file.info.table$input.files[tmp.lines[i]]))
        if (is.list(w)) w <- convertListTocompData(w)
        if (!is.null(w)) tmp.uID <- c(tmp.uID, info.parameters(w)$uID)
      }
    }
    if (length(unique(tmp.uID)) > 1) {
      stop("Different unique identifiers for data sets with the same name, number of 
           samples and replicate ID.")
    }
  }
  
  ## Then check that all input files with the same dataset name (i.e., all input files) 
  ## only differ in the nbr samples and repl.id
  nex.files <- NULL
  ii <- 1
  while (!file.exists(as.character(file.info.table$input.files[ii]))) {
    nex.files <- c(nex.files, ii)
    ii <- ii + 1
  }
  w <- readRDS(as.character(file.info.table$input.files[ii]))
  if (is.list(w)) w <- convertListTocompData(w)
  tmp.ref <- info.parameters(w)
  tmp.ref[c("uID", "samples.per.cond", "repl.id")] <- NULL
  
  if (ii != length(file.info.table$input.files)) {
    for (i in (ii + 1):length(file.info.table$input.files)) {
      if (file.exists(as.character(file.info.table$input.files[i]))) {
        w <- readRDS(as.character(file.info.table$input.files[i]))
        if (is.list(w)) w <- convertListTocompData(w)
        tmp.comp <- info.parameters(w)
        tmp.comp[c("uID", "samples.per.cond", "repl.id")] <- NULL
        if (!(all(tmp.ref %in% tmp.comp) & all(tmp.comp %in% tmp.ref))) {
          print(tmp.ref[!(tmp.ref %in% tmp.comp)])
          print(tmp.comp[!(tmp.comp %in% tmp.ref)])
          stop("The same data set name is used to refer to data sets with different 
           simulation parameters")
        }
      } else {
        nex.files <- c(nex.files, i)
      }
    }
  }
}

## Comparison functions
# Assign plot colors to a collection of differential expression methods
# 
# Assign a plot color to each differential expression method present in the input vector. This color will be used for all graphical representations in the method comparison. The colors are taken from a pre-defined collection. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparison}} function.
# 
# @param methods Vector of methods that will be assigned plot colors.
# @return List with the assigned plot colors for all methods present in \code{methods}.
# @author Charlotte Soneson
definePlotColors <- function(methods) {
  available.plot.colors = c('red', 'cyan', 'blue', 'green', 'orange', 'violet', 'brown',
                            'deeppink', 'forestgreen', 'aquamarine4', 'cadetblue3', 'grey', 
                            'bisque', 'darkkhaki', 'darkgoldenrod3', 'pink', 'plum',
                            'tomato', 'wheat', 'azure4', 'blue4', 'blueviolet', 'brown4', 
                            'burlywood3', 'burlywood4', 'chartreuse', 'chocolate',
                            'cornflowerblue', 'cornsilk3', 'darkolivegreen2', 'darkorange2', 
                            'darksalmon', 'darkseagreen', 'dodgerblue', 'firebrick1',
                            'lavender', 'lemonchiffon2', 'hotpink', 'lightskyblue2', 
                            'mediumaquamarine', 'olivedrab2', 'palegreen')
  plot.colors <- list()
  for (i in seq_len(length(methods))) {
    plot.colors[sort(methods)[i]] <- available.plot.colors[i]
  }
  return(plot.colors)
}

# Plot ROC curves
# 
# Plot receiver operating characteristic (ROC) curves, depicting the relationship between the false positive rate (x-axis) and the true positive rate (y-axis) as the significance threshold for differential expression is changed. The genes are ordered by the \code{score} column of the \code{result.table} data frame available in the result objects returned by \code{\link{runDiffExp}}. A well-performing differential expression method has a sharply increasing ROC curve, with a large area under the curve (AUC). This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparison}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which ROC curves will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which ROC curves will be constructed
# @author Charlotte Soneson
makeROCcurves = function(setup.parameters, sel.nbrsamples, sel.repl) {
  ## Compute and plot ROC curves for the selected data set(s).
  ## One figure for each sel.nbrsamples, one panel for each sel.repl
  nbr.cols <- 3
  nbr.rows <- ceiling((length(sel.repl) + 1)/nbr.cols)
  
  for (i in seq_len(length(sel.nbrsamples))) {
    incl.legend <- FALSE
    graphics::par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
                  cex.axis = 1.5, las = 1, cex.main = 1.75, mar = c(5, 5, 4, 2))
    all.plot.methods <- NULL
    all.plot.colors <- NULL
    for (j in seq_len(length(sel.repl))) {
      plot.colors <- NULL
      plot.methods <- NULL
      tmp.k <- 1
      for (k in seq_len(length(setup.parameters$incl.de.methods))) {
        ## Extract results from selected DE method, nbr samples and replicate
        idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
                            setup.parameters$incl.de.methods[k])
        idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
                            sel.nbrsamples[i])
        idx3 <- findFileIdx(setup.parameters$file.info, "repl", 
                            sel.repl[j])
        idx <- intersect(intersect(idx1, idx2), idx3)
        if (length(idx) != 0) {
          ## Read data
          X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
          if (is.list(X)) X <- convertListTocompData(X)
          the.scores <- result.table(X)$score
          if (!all(variable.annotations(X)$differential.expression == 0)) {
            incl.legend <- TRUE
            idx_not_na <- which(!is.na(the.scores))
            the.prediction <- prediction(the.scores[idx_not_na],
                                         variable.annotations(X)$differential.expression[idx_not_na])
            the.performance <- performance(the.prediction,
                                           measure = 'tpr',
                                           x.measure = 'fpr')
            fpr <- unlist(the.performance@x.values)
            tpr <- unlist(the.performance@y.values)
            
            plot.colors <- c(plot.colors,
                            setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]])
            plot.methods <- c(plot.methods, setup.parameters$incl.de.methods[k])
            
            ## Create vectors for legend (combine info across all replicates)
            if (length(which(all.plot.methods == setup.parameters$incl.de.methods[k])) == 0) {
              all.plot.methods <- c(all.plot.methods, setup.parameters$incl.de.methods[k])
              all.plot.colors <- c(all.plot.colors, setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]])
            }
            
            if (tmp.k == 1) {
              if (is.na(sel.repl[j])) {
                if (is.na(sel.nbrsamples[i])) {
                  maintext <- ""
                } else {
                  maintext <- paste(sel.nbrsamples[i], "samples/condition")
                }
              } else {
                if (is.na(sel.nbrsamples[i])) {
                  maintext <- paste("Replicate", sel.repl[j])
                } else {
                  maintext <- paste("Replicate", sel.repl[j], ", ", sel.nbrsamples[i], 
                                    "samples/condition")
                }
              }
              plot(fpr, tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1),
                   col = setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]],
                   main = maintext, 
                   lwd = 1.5, 
                   xlab = 'False positive rate',
                   ylab = 'True positive rate')
            } else {
              graphics::lines(
                fpr, tpr, 
                col = setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]], 
                lwd = 1.5)
            }
          } 
          tmp.k <- tmp.k + 1
        }
      }
    }
    if (incl.legend) {
      plot(0:1, 0:1, type = 'n', axes = FALSE, xlab = '', ylab = '')
      graphics::legend('topleft', legend = all.plot.methods, 
                       col = all.plot.colors, lty = 1, lwd = 2, 
                       cex = ifelse(length(all.plot.methods) < 24, 1.5, 
                                    1.3/length(all.plot.methods)*24))
    } else {
      cat("No truly differentially expressed genes found. ")
    }
  }
}

# Plot false discovery curves
# 
# Plot false discovery curves, depicting the number of false discoveries encountered as the significance threshold for differential expression is changed. The genes are ordered by the \code{score} column of the \code{result.table} data frame available in the result objects returned by \code{\link{runDiffExp}}. A well-performing differential expression method has a slowly increasing false discovery curve, meaning that as more and more genes are considered differentially expressed, only few of them are false discoveries. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparison}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the false discovery curves will be constructed. 
# @param sel.repl The replicate numbers (the instances of a given simulation setting) for which the false discovery curves will be constructed.
# @author Charlotte Soneson
makeFalseDiscoveryCurves <- function(setup.parameters, sel.nbrsamples, sel.repl, legend.outside = TRUE) {
  nbr.cols <- 3
  nbr.rows <- ceiling((length(sel.repl) + 1)/nbr.cols)
  
  ##sel.repl <- sort(sel.repl)
  
  for (i in seq_len(length(sel.nbrsamples))) {
    incl.legend <- FALSE
    graphics::par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
                  cex.axis = 1.5, las = 1, cex.main = 1.75, mar = c(5, 5, 4, 2), 
                  mgp = c(3.5, 1, 0))
    all.plot.colors <- NULL
    all.plot.methods <- NULL
    for (j in seq_len(length(sel.repl))) {
      plot.colors <- NULL
      plot.methods <- NULL
      tmp.k <- 1
      for (k in seq_len(length(setup.parameters$incl.de.methods))) {
        idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
                            setup.parameters$incl.de.methods[k])
        idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
                            sel.nbrsamples[i])
        idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl[j])
        idx <- intersect(intersect(idx1, idx2), idx3)
        if (length(idx) != 0) {
          X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
          if (is.list(X)) X <- convertListTocompData(X)
          the.scores <- result.table(X)$score
          if (!all(variable.annotations(X)$differential.expression == 0)) {
            incl.legend <- TRUE
            the.scores.ordered <- the.scores[order(the.scores,
                                                  decreasing = TRUE)]
            trulyDE.ordered <- variable.annotations(X)$differential.expression[order(the.scores,
                                                                                   decreasing = TRUE)]
            the.scores.unique <- sort(unique(the.scores), decreasing = TRUE)
            xcoord <- cumsum(as.vector(table(-the.scores.ordered)))
            ycoord <- cumsum(1 - trulyDE.ordered)[cumsum(as.vector(table(-the.scores.ordered)))]
            
            plot.colors <- c(plot.colors,
                             setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]])
            plot.methods <- c(plot.methods, setup.parameters$incl.de.methods[k])
            
            ## Create vectors for legend (combine info across all replicates)
            if (length(which(all.plot.methods == setup.parameters$incl.de.methods[k])) == 0) {
              all.plot.methods <- c(all.plot.methods, setup.parameters$incl.de.methods[k])
              all.plot.colors <- c(all.plot.colors, setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]])
            }
            
            if (all(xcoord < setup.parameters$fdc.maxvar)) {
              tmp.fdc.maxvar <- max(xcoord)
            } else {
              tmp.fdc.maxvar <- setup.parameters$fdc.maxvar
            }
            
            if (tmp.k == 1) {
              if (is.na(sel.repl[j])) {
                if (is.na(sel.nbrsamples[i])) {
                  maintext <- ""
                } else {
                  maintext <- paste(sel.nbrsamples[i], "samples/condition")
                }
              } else {
                if (is.na(sel.nbrsamples[i])) {
                  maintext <- paste("Replicate", sel.repl[j])
                } else {
                  maintext <- paste("Replicate", sel.repl[j], ", ", sel.nbrsamples[i], 
                                    "samples/condition")
                }
              }
              
              suppressWarnings({plot(xcoord[seq_len(min(which(xcoord >= tmp.fdc.maxvar)))],
                                     ycoord[seq_len(min(which(xcoord >= tmp.fdc.maxvar)))],
                                     type = 'l', xlim = c(0, tmp.fdc.maxvar), log = 'y', las = 1,
                                     lwd = 1.5,  #1.5
                                     ylim = c(1, tmp.fdc.maxvar),
                                     col = setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]], 
                                     xlab = 'Number of selected genes',
                                     ylab = 'Number of false discoveries',
                                     main = maintext)})
            } else {
              suppressWarnings({
                graphics::lines(xcoord[seq_len(min(which(xcoord >= tmp.fdc.maxvar)))],
                                ycoord[seq_len(min(which(xcoord >= tmp.fdc.maxvar)))],
                                col = setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]], 
                                lwd = 1.5)})   #1.5
            }
          } 
          tmp.k <- tmp.k + 1
        }
      }
    }
    if (incl.legend) {
      if (legend.outside) {
        plot(0:1, 0:1, type = 'n', axes = FALSE, xlab = '', ylab = '')
        graphics::legend('topleft', legend = all.plot.methods, 
                         col = all.plot.colors, lty = 1, lwd = 2, 
                         cex = ifelse(length(all.plot.methods) < 24, 1.5, 
                                      1.3/length(all.plot.methods)*24))
      } else {
        graphics::legend('bottomright', legend = all.plot.methods, 
                         col = all.plot.colors, lty = 1, lwd = 2, 
                         cex = ifelse(length(all.plot.methods) < 24, 1.5, 
                                      1.3/length(all.plot.methods)*24))
      }
    } else {
      cat("No truly differentially expressed genes found.")
    }
  }
}

# Construct MA plots
# 
# Construct MA plots, depicting the log-fold change (M) vs the average expression level (A) for all genes and indicating the genes that are determined to be differentially expressed. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the MA plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the MA plots will be constructed
# @author Charlotte Soneson
plotMASignificant <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  nbr.cols <- 3
  nbr.rows <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols)
  for (i in seq_len(length(sel.nbrsamples))) {
    graphics::par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
                  cex.axis = 1.5, las = 1, cex.main = 1.2, 
                  cex.main = 1.75, mar = c(5, 5, 4, 2))
    for (k in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods",
                          setup.parameters$incl.de.methods[k])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples",
                          sel.nbrsamples[i])
      idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
      idx <- intersect(intersect(idx1, idx2), idx3)
      if (length(idx) != 0) {
        X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
        if (is.list(X)) X <- convertListTocompData(X)
        if (is.null(variable.annotations(X)$A.value) | is.null(variable.annotations(X)$M.value)) {
          A.value <- computeAval(count.matrix(X), sample.annotations(X)$condition)
          M.value <- computeMval(count.matrix(X), sample.annotations(X)$condition)
        } else {
          A.value <- variable.annotations(X)$A.value
          M.value <- variable.annotations(X)$M.value
        }
        if ('adjpvalue' %in% colnames(result.table(X))) {
          significant <- which(result.table(X)$adjpvalue < setup.parameters$ma.threshold)
        } else if ('FDR' %in% colnames(result.table(X))) {
            significant <- which(result.table(X)$FDR < setup.parameters$ma.threshold)
        } else {
          significant <- NULL
        }
        maplotcols <- rep('black', length(A.value))
        maplotcols[significant] <- setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]]
        if (!is.null(significant)) {
          plot(A.value, M.value, cex = 1, pch = 20, main = setup.parameters$incl.de.methods[k], 
               xlab = "Average expression", 
               ylab = "Log fold change", 
               col = maplotcols)
        } else {
          plot(1, type = "n", xlab = "Average expression", 
               ylab = "Log fold change",
               main = setup.parameters$incl.de.methods[k], 
               xaxt = "n", yaxt = "n")
        }
      }
    }
  }
}

# Plot the gene score vs the average expression level of the gene
# 
# Plot the gene score vs the average expression level (A) for all genes. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the MA plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the plots will be constructed
# @author Charlotte Soneson
plotScoreVsExpr <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  nbr.cols <- 3
  nbr.rows <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols)
  for (i in seq_len(length(sel.nbrsamples))) {
    graphics::par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
                  cex.axis = 1.5, las = 1, cex.main = 1.2, 
                  cex.main = 1.75, mar = c(5, 5, 4, 2))
    for (k in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
                          setup.parameters$incl.de.methods[k])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
                          sel.nbrsamples[i])
      idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
      idx <- intersect(intersect(idx1, idx2), idx3)
      if (length(idx) != 0) {
        X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
        if (is.list(X)) X <- convertListTocompData(X)
        if (is.null(variable.annotations(X)$A.value)) {
          A.value <- computeAval(count.matrix(X), sample.annotations(X)$condition)
        } else {
          A.value <- variable.annotations(X)$A.value
        }
        score <- result.table(X)$score
        linecol <- setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]]
        plot(A.value, score, cex = 0.75, pch = 20, main = setup.parameters$incl.de.methods[k], 
             xlab = "Average expression", 
             ylab = "Score")
        loessline <- stats::loess(score ~ A.value)
        xval <- seq(min(A.value), max(A.value), length.out = 500)
        graphics::lines(xval, stats::predict(loessline, xval), col = linecol, lwd = 5) 
      }
    }
  }
}

# Plot the gene score vs an 'outlier score'
# 
# Plot the gene score vs an outlier score, defined as the ratio between A and B, where A is the difference between the highest and the second highest normalized count, and B is the difference between the highest and the lowest normalized count. The outlier score is averaged across all samples. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the MA plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the MA plots will be constructed
# @author Charlotte Soneson
# plotScoreVsOutlierEvidence <- function(setup.parameters, sel.nbrsamples, sel.repl) {
#   nbr.cols <- 3
#   nbr.rows <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols)
#   for (i in seq_len(length(sel.nbrsamples))) {
#     par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
#         cex.axis = 1.5, las = 1, cex.main = 1.2, cex.main = 1.75, mar = c(5, 5, 4, 2))
#     for (k in seq_len(length(setup.parameters$incl.de.methods))) {
#       idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
#                           setup.parameters$incl.de.methods[k])
#       idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
#                           sel.nbrsamples[i])
#       idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
#       idx <- intersect(intersect(idx1, idx2), idx3)
#       if (length(idx) != 0) {
#         X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
#         if (is.list(X)) X <- convertListTocompData(X)
#         ## Normalize counts (compute pseudocounts)
#         nf <- calcNormFactors(count.matrix(X))
#         norm.factors <- nf * colSums(count.matrix(X))
#         common.libsize <- exp(mean(log(colSums(count.matrix(X)))))
#         pseudocounts <- sweep(count.matrix(X) + 0.5, 2, norm.factors, '/') * common.libsize
#         ## Compute the outlier score
#         tmp <- t(apply(count.matrix(X)[, sample.annotations(X)$condition == 
#                                          levels(factor(sample.annotations(X)$condition))[1]], 
#                        1, sort))
#         outlier.score.1.1 <- (tmp[, ncol(tmp)] - tmp[, (ncol(tmp) - 1)])/(tmp[, ncol(tmp)] - 
#                                                                             tmp[, 1])
#         outlier.score.1.2 <- (tmp[, 2] - tmp[, 1])/(tmp[, ncol(tmp)] - tmp[, 1])
#         outlier.score.1 <- sqrt(outlier.score.1.1*outlier.score.1.2)
#         outlier.score.1[is.na(outlier.score.1)] <- 0
#         
#         tmp <- t(apply(count.matrix(X)[, sample.annotations(X)$condition == 
#                                          levels(factor(sample.annotations(X)$condition))[2]], 
#                        1, sort))
#         outlier.score.2.1 <- (tmp[, ncol(tmp)] - tmp[, (ncol(tmp) - 1)])/(tmp[, ncol(tmp)] - 
#                                                                           tmp[, 1])
#         outlier.score.2.2 <- (tmp[, 2] - tmp[, 1])/(tmp[, ncol(tmp)] - tmp[, 1])
#         outlier.score.2 <- sqrt(outlier.score.2.1*outlier.score.2.2)
#         outlier.score.2[is.na(outlier.score.2)] <- 0
#         outlier.score <- 0.5*(outlier.score.1 + outlier.score.2)
#         
#         ## Count number of outliers
#         if (length(variable.annotations(X)) != 0) {
#           total.outliers <- apply(variable.annotations(X)[match(rownames(result.table(X)), 
#                                                                 rownames(variable.annotations(X))), 
#                                                           which(colnames(variable.annotations(X)) %in% 
#                                                                   c("n.random.outliers.up.S1", 
#                                                                     "n.random.outliers.up.S2", 
#                                                                     "n.random.outliers.down.S1", 
#                                                                     "n.random.outliers.down.S2", 
#                                                                     "n.single.outliers.up.S1", 
#                                                                     "n.single.outliers.up.S2", 
#                                                                     "n.single.outliers.down.S1", 
#                                                                     "n.single.outliers.down.S2"))], 1, sum)
#         } else {
#           total.outliers <- rep(0, length(score))
#         }
#         
# #        boxplot(outlier.score ~ total.outliers)
#         score <- result.table(X)$score
#         linecol <- setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]]
#         plot(outlier.score, score, cex = 0.75, pch = 20, 
#              main = setup.parameters$incl.de.methods[k], 
#              xlab = "Outlier score", 
#              ylab = "Score")
#         loessline <- loess(score ~ outlier.score)
#         xval <- seq(min(outlier.score), max(outlier.score), length.out = 500)
#         lines(xval, predict(loessline, xval), col = linecol, lwd = 5) 
#       }
#     }
#   }
# }

findGenesAllZero <- function(count.matrix, conditions) {
  ## Find genes that are zero for all samples in one condition
  cm1 <- count.matrix[, which(conditions == levels(factor(conditions))[1])]
  cm2 <- count.matrix[, which(conditions == levels(factor(conditions))[2])]
  zero1 <- which(apply(cm1, 1, sum) == 0)
  zero2 <- which(apply(cm2, 1, sum) == 0)
  return(list(zero1 = zero1, zero2 = zero2))
}

computeSignal <- function(count.matrix, conditions, signal.measure) {
  ## Normalize counts (compute pseudocounts) and log-transform
  nf <- calcNormFactors(count.matrix)
  norm.factors <- nf * colSums(count.matrix)
  common.libsize <- exp(mean(log(colSums(count.matrix))))
  pseudocounts <- sweep(count.matrix + 0.5, 2, norm.factors, '/') * common.libsize
  log2.pseudocounts <- log2(pseudocounts)
  
  ## If the second line below is uncommented, compute the "log2 pseudocounts" for 
  ## a matrix with all zeros (the same normalization factors as above) and 
  ## subtract from the actual log2 pseudocounts
  tmpmat <- matrix(0, nrow(count.matrix), ncol(count.matrix))
  ##tmpmat <- log2(sweep(tmpmat + 0.5, 2, norm.factors, "/") * common.libsize)
  
  ## Compute difference between the values
  log2.pseudocounts.diff <- log2.pseudocounts - tmpmat
  
  if (signal.measure == "snr") {
    ## Split into two matrices depending on condition. Average within each condition for each gene
    signal1 <- 
      apply(log2.pseudocounts.diff[, which(conditions == 
                                             levels(factor(conditions))[1])], 
            1, mean)/apply(log2.pseudocounts.diff[, which(conditions == 
                                                            levels(factor(conditions))[1])], 
                           1, stats::sd)
    signal2 <- 
      apply(log2.pseudocounts.diff[, which(conditions == 
                                             levels(factor(conditions))[2])], 
            1, mean)/apply(log2.pseudocounts.diff[, which(conditions == 
                                                            levels(factor(conditions))[2])], 
                           1, stats::sd)
  } else if (signal.measure == "mean") {
    signal1 <- apply(log2.pseudocounts.diff[, which(conditions == 
                                                      levels(factor(conditions))[1])], 
                  1, mean)
    signal2 <- apply(log2.pseudocounts.diff[, which(conditions == 
                                                      levels(factor(conditions))[2])], 
                  1, mean)
  }
  
  return(list(signal1 = signal1, signal2 = signal2))
}

# Plot the gene score vs the "signal", for genes expressed in only one condition
# 
# For genes expressed in only one condition, plot the gene score vs the average expression level in that condition. All replicates are merged in one plot. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the plots will be constructed
# @author Charlotte Soneson
plotSignalForZeroCounts <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  theplotcols <- c("black", "red", "green", "blue", "grey", "cyan", "forestgreen",
                   "orange", "brown", "violet", "deeppink", "aquamarine4", 
                   "cadetblue3", "bisque", "darkkhaki", "darkgoldenrod3", "pink",
                   "plum", "tomato", "yellow", "azure4", "blue4", "blueviolet", 
                   "chartreuse", "chocolate", "cornflowerblue", "darkorange2", 
                   "darksalmon", "lemonchiffon2", "palegreen", "lightskyblue2")
  
  nbr.cols <- 3
  nbr.rows <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols)
  for (i in seq_len(length(sel.nbrsamples))) {
    graphics::par(mfrow = c(nbr.rows, nbr.cols), cex.lab = 2, 
                  cex.axis = 1.5, las = 1, cex.main = 1.2, 
                  cex.main = 1.75, mar = c(5, 5, 4, 2))
    for (k in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
                          setup.parameters$incl.de.methods[k])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
                          sel.nbrsamples[i])
      tmpx <- c()
      tmpy <- c()
      plotcol <- c()
      for (j in seq_len(length(sel.repl))) {
        idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl[j])
        idx <- intersect(intersect(idx1, idx2), idx3)
        if (length(idx) != 0) {
          X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
          if (is.list(X)) X <- convertListTocompData(X)
          zerocounts <- findGenesAllZero(count.matrix(X), sample.annotations(X)$condition)
          signals <- computeSignal(count.matrix(X), sample.annotations(X)$condition, 
                                   signal.measure = setup.parameters$signal.measure)
          tmpx <- c(tmpx, signals$signal1[zerocounts$zero2], signals$signal2[zerocounts$zero1])
          tmpy <- c(tmpy, result.table(X)$score[zerocounts$zero2], 
                    result.table(X)$score[zerocounts$zero1])
          plotcol <- c(plotcol, rep(theplotcols[j], 
                                    length(zerocounts$zero1) + length(zerocounts$zero2)))
        }
      }
      if ((length(tmpx) != 0) && (!(all(is.na(tmpy))))) {
        plot(tmpx, tmpy, xlab = setup.parameters$signal.measure, 
             ylab = "Score", cex = 1, pch = 20, col = plotcol, 
             main = setup.parameters$incl.de.methods[k])
        ##linecol <- setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]]
        ##loessline <- loess(tmpy ~ tmpx, span = 0.25)
        ##xval <- seq(min(tmpx), max(tmpx), length.out = 500)
        ##lines(xval, predict(loessline, xval), col = linecol, lwd = 5) 
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
             main = setup.parameters$incl.de.methods[k])
      }
    }
  }
}

# computeM <- function(count.matrix, conditions) {
#   nf <- calcNormFactors(count.matrix)
#   norm.factors <- nf * colSums(count.matrix)
#   common.libsize <- exp(mean(log(colSums(count.matrix))))
#   pseudocounts <- sweep(count.matrix + 0.5, 2, norm.factors, '/') * common.libsize
#   log2.pseudocounts <- log2(pseudocounts)
#   M.value <- apply(log2.pseudocounts[, which(conditions == levels(factor(conditions))[2])], 
#                    1, mean) - 
#     apply(log2.pseudocounts[, which(conditions == levels(factor(conditions))[1])], 
#           1, mean)
#   return(M.value)
# }

# computeA <- function(count.matrix, conditions) {
#   nf <- calcNormFactors(count.matrix)
#   norm.factors <- nf * colSums(count.matrix)
#   common.libsize <- exp(mean(log(colSums(count.matrix))))
#   pseudocounts <- sweep(count.matrix + 0.5, 2, norm.factors, '/') * common.libsize
#   log2.pseudocounts <- log2(pseudocounts)
#   A.value <- 0.5*(apply(log2.pseudocounts[, which(conditions == levels(factor(conditions))[2])], 
#                         1, mean) + 
#                     apply(log2.pseudocounts[, which(conditions == 
#                                                       levels(factor(conditions))[1])], 
#                           1, mean))
#   return(A.value)
# }

computeFDR <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues < signthreshold), 
                   which(trueDElabels == 0)))/length(which(adjpvalues < signthreshold))
}

computeTPR <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues < signthreshold),
                   which(trueDElabels == 1)))/length(which(trueDElabels == 1))
}

computeTP <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues < signthreshold),
                   which(trueDElabels == 1)))
}

computeFP <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues < signthreshold),
                   which(trueDElabels == 0)))
}

computeTN <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues >= signthreshold),
                   which(trueDElabels == 0)))
}

computeFN <- function(adjpvalues, trueDElabels, signthreshold) {
  length(intersect(which(adjpvalues >= signthreshold),
                   which(trueDElabels == 1)))
}

computeMCC <- function(adjpvalues, trueDElabels, signthreshold) {
  tp <- length(intersect(which(adjpvalues < signthreshold),
                          which(trueDElabels == 1)))
  fp <- length(intersect(which(adjpvalues < signthreshold),
                          which(trueDElabels == 0)))
  tn <- length(intersect(which(adjpvalues >= signthreshold),
                          which(trueDElabels == 0)))
  fn <- length(intersect(which(adjpvalues >= signthreshold),
                          which(trueDElabels == 1)))
  n <- tp + fp + tn + fn
  s <- (tp + fn)/n
  p <- (tp + fp)/n
  (tp/n - s*p)/sqrt(p*s*(1 - s)*(1 - p))
}

computeTypeIerror <- function(pvalues, trueDElabels, signthreshold) {
  length(intersect(which(pvalues < signthreshold),
                   which(trueDElabels == 0)))/length(which(trueDElabels == 0))
  
}

findFileIdx <- function(file.table, column, value) {
  if (is.na(value)) {
    return(which(is.na(file.table[, column])))
  } else {
    return(which(file.table[, column] == value))
  }
}

# Plot the observed FDR vs the average expression level of the gene
# 
# Plot the observed FDR vs the average expression level (A) for all genes. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the MA plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the MA plots will be constructed
# @author Charlotte Soneson
plotFDRVsExpr <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  for (i in seq_len(length(sel.nbrsamples))) {
    all.methods <- NULL
    obs.fdr.total <- NULL
    doplot <- FALSE
    idx1 <- findFileIdx(setup.parameters$file.info, "nbr.samples", sel.nbrsamples[i])
    for (k in seq_len(length(setup.parameters$incl.de.methods))) {
      obs.fdr.method <- NULL
      idx2 <- intersect(idx1, findFileIdx(setup.parameters$file.info, "de.methods", 
                                          setup.parameters$incl.de.methods[k]))
      for (j in seq_len(length(sel.repl))) {
        idx <- intersect(idx2, findFileIdx(setup.parameters$file.info, "repl", 
                                           sel.repl[j]))
        if (length(idx) != 0) {
          X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
          if (is.list(X)) X <- convertListTocompData(X)
          if (!all(variable.annotations(X)$differential.expression == 0)) {
            if (is.null(variable.annotations(X)$A.value)) {
              A.value <- computeAval(count.matrix(X), sample.annotations(X)$condition)
            } else {
              A.value <- variable.annotations(X)$A.value
            }
            if ('adjpvalue' %in% colnames(result.table(X))) {
              adjp <- result.table(X)$adjpvalue
            } else if ('FDR' %in% colnames(result.table(X))) {
              adjp <- result.table(X)$FDR
            } else {
              adjp <- NULL
            }
            if (!is.null(adjp)) {
              if (!(setup.parameters$incl.de.methods[k] %in% all.methods)) {
                all.methods <- c(all.methods, setup.parameters$incl.de.methods[k])
              }
              doplot <- TRUE
              tmp <- data.frame(status = variable.annotations(X)$differential.expression, 
                                adjp = adjp)
              ## Split the genes by expression level 
              r <- range(A.value)
              expr.bin <- cut(A.value, seq(r[1], r[2], length.out = 11), 
                              include.lowest = TRUE)
              levels(expr.bin) <- c("(0,10]", "(10,20]", "(20,30]", "(30,40]", 
                                    "(40,50]", "(50,60]", "(60,70]", "(70,80]", 
                                    "(80,90]", "(90,100]")
              ## Split the fdr vector and compute the observed FDR
              split.tmp <- split(tmp, expr.bin)
              obs.fdr <- vapply(split.tmp, function(x, thr = setup.parameters$fdr.threshold) {
                computeFDR(x$adjp, x$status, thr)
              }, FUN.VALUE = NA_real_)
              obs.fdr.method <- rbind(obs.fdr.method, data.frame(bin = names(obs.fdr), 
                                                                 fdr = obs.fdr, 
                                                                 method = rep(setup.parameters$incl.de.methods[k], length(obs.fdr))))
            }
          }
        }
      }
      if ((!is.null(obs.fdr.method)) && (!all(is.na(obs.fdr.method$fdr)))) {
        obs.fdr.total <- rbind(obs.fdr.total, obs.fdr.method)
      }
    }
    #}
    if (doplot && is.data.frame(obs.fdr.total)) {
      obs.fdr.total$method <- factor(obs.fdr.total$method, 
                                     levels = setup.parameters$incl.de.methods)
      
      suppressWarnings({
        print(ggplot(obs.fdr.total, aes_string(x = "bin", y = "fdr", fill = "method")) + 
                geom_boxplot(outlier.size = 0) +
                scale_fill_manual(values = unlist(setup.parameters$specified.colors), 
                                  name = "", guide = FALSE) + 
                facet_wrap(~ method, ncol = 3, drop = FALSE, scales = "free_x") + 
                geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 1.5) + 
                xlab("Bin") + 
                ylab("FDR") + 
                theme_bw() + 
                scale_y_continuous(limits = c(-0, 1)) + 
                theme(strip.text.x = element_text(size = 15), 
                      axis.text.y = element_text(size = 15),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), 
                      axis.title.x = element_text(size = 15),
                      axis.title.y = element_text(size = 15)))
      })
      
    } else {
      cat ("No truly differentially expressed genes found, or no differentially expressed genes found by any method. ")
    }
  }
}

# Plot distribution of gene scores as a function of the number of outliers
# 
# Plot the distribution of the scores assigned to the genes, grouped by the total number of outlier counts introduced for the genes. 
# 
# @param setup.parameters List of parameters (internal).
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the enrichment plots will be constructed
# @param sel.repl The replicate numbers (instances of a given simulation setting) for which the enrichment plots will be constructed
# @author Charlotte Soneson
plotScoreVsOutliers <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  nbr.cols <- 3
  nbr.rows <- ceiling(length(setup.parameters$incl.de.methods)/nbr.cols)
  for (i in seq_len(length(sel.nbrsamples))) {
    graphics::par(mfrow = c(nbr.rows, nbr.cols), mar = c(9, 6, 4, 2), mgp = c(4, 1, 0), 
                  cex.axis = 1.5, las = 1, cex.main = 1.75, cex.lab = 2)
    for (k in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods", 
                          setup.parameters$incl.de.methods[k])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples", 
                          sel.nbrsamples[i])
      idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
      idx <- intersect(intersect(idx1, idx2), idx3)
      if (length(idx) != 0) {
        X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
        if (is.list(X)) X <- convertListTocompData(X)
        score <- result.table(X)$score
        if (length(variable.annotations(X)) != 0) {
          total.outliers <- apply(variable.annotations(X)[match(rownames(result.table(X)), 
                                                                rownames(variable.annotations(X))), 
                                                          which(colnames(variable.annotations(X)) %in% 
                                                                  c("n.random.outliers.up.S1", 
                                                                    "n.random.outliers.up.S2", 
                                                                    "n.random.outliers.down.S1", 
                                                                    "n.random.outliers.down.S2", 
                                                                    "n.single.outliers.up.S1", 
                                                                    "n.single.outliers.up.S2", 
                                                                    "n.single.outliers.down.S1", 
                                                                    "n.single.outliers.down.S2"))], 1, sum)
        } else {
          total.outliers <- rep(0, length(score))
        }
        
        ## Violin plots
        keep <- which(total.outliers %in% 
                        names(table(total.outliers))[which(table(total.outliers) > 3)])
        total.outliers <- total.outliers[keep]
        score <- score[keep]
                
        w <- split(score, total.outliers)
        names(w)[1] <- "x"
        
        do.call(vioplot, c(lapply(w, stats::na.omit), 
                           list(horizontal = FALSE, names = names(table(total.outliers))),
                           col = setup.parameters$specified.colors[[setup.parameters$incl.de.methods[k]]]))
        graphics::title(main = setup.parameters$incl.de.methods[k], 
                        xlab = "Total number of outliers", 
                        ylab = "Score")
                
        graphics::axis(1, at = seq_len(length(unique(total.outliers))), 
                       labels = paste("n =", table(total.outliers)[match(levels(factor(total.outliers)), 
                                                                         names(table(total.outliers)))]), 
                       tick = FALSE, line = 1, lwd = 0)
        
      }
    }
  }
}

# Construct a result table with selected performance measures for later plotting
# 
# Construct a result table containing one or several of the type I error, the fraction of significant genes, the false discovery rate (FDR), the true positive rate (TPR), Matthew's correlation coefficient (MCC), and the area under the ROC curve (AUC). This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# The result table returned from this function may contain "missing" lines, for example, there may be some differential expression methods that have not been applied to all data sets. This can be fixed by passing the result table through the \code{\link{padResultTable}} function, which will fill in all such missing lines with NA values.
# 
# @param setup.parameters List of parameters (internal). The \code{comparisons} entry defines the performance measures that will be included in the result table.
# @return A result table, containing the selected performance measures for each differential expression method and for each of the included data sets (the collection of data sets to include is defined by the \code{setup.parameters} list). 
# @author Charlotte Soneson
createResultTable <- function(setup.parameters) {
  if ('typeIerror' %in% setup.parameters$comparisons) {tIe.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('fracsign' %in% setup.parameters$comparisons) {fsn.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('nbrsign' %in% setup.parameters$comparisons) {nsn.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('nbrtpfp' %in% setup.parameters$comparisons) {
    tp.vec <- rep(NA, nrow(setup.parameters$file.info))
    fp.vec <- rep(NA, nrow(setup.parameters$file.info))
    tn.vec <- rep(NA, nrow(setup.parameters$file.info))
    fn.vec <- rep(NA, nrow(setup.parameters$file.info))
  }
  if ('fdr' %in% setup.parameters$comparisons) {fdr.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('tpr' %in% setup.parameters$comparisons) {tpr.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('auc' %in% setup.parameters$comparisons) {auc.vec <- rep(NA, nrow(setup.parameters$file.info))}
  if ('mcc' %in% setup.parameters$comparisons) {mcc.vec <- rep(NA, nrow(setup.parameters$file.info))}
  
  for (i in seq_len(nrow(setup.parameters$file.info))) {
    X <- readRDS(as.character(setup.parameters$file.info$input.files[i]))
    if (is.list(X)) X <- convertListTocompData(X)
    if ('typeIerror' %in% setup.parameters$comparisons && 'pvalue' %in% colnames(result.table(X))) {
      tIe.vec[i] <- computeTypeIerror(result.table(X)$pvalue, 
                                      variable.annotations(X)$differential.expression,
                                      setup.parameters$typeI.threshold)
    }
    if (!all(variable.annotations(X)$differential.expression == 0)) {
      if ('auc' %in% setup.parameters$comparisons) {
        the.scores <- result.table(X)$score
        idx_not_na <- which(!is.na(the.scores))
        the.prediction <- prediction(the.scores[idx_not_na], 
                                     variable.annotations(X)$differential.expression[idx_not_na])
        auc.vec[i] <- unlist(performance(the.prediction, measure = 'auc')@y.values)
      }
      
      if ('mcc' %in% setup.parameters$comparisons) {
        if ('FDR' %in% colnames(result.table(X))) {
          mcc.vec[i] <- computeMCC(result.table(X)$FDR,
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$mcc.threshold)
        } else {
          if ('adjpvalue' %in% colnames(result.table(X))) {
            mcc.vec[i] <- computeMCC(result.table(X)$adjpvalue,
                                     variable.annotations(X)$differential.expression,
                                     setup.parameters$mcc.threshold)
          }
        }
      }
      
      if ('fdr' %in% setup.parameters$comparisons) {
        if ('FDR' %in% colnames(result.table(X))) {
          fdr.vec[i] <- computeFDR(result.table(X)$FDR, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$fdr.threshold)
        } else {
          if ('adjpvalue' %in% colnames(result.table(X))) {
            fdr.vec[i] <- computeFDR(result.table(X)$adjpvalue, 
                                     variable.annotations(X)$differential.expression,
                                     setup.parameters$fdr.threshold)
          }
        }
      }
      
      if ('tpr' %in% setup.parameters$comparisons) {
        if ('FDR' %in% colnames(result.table(X))) {
          tpr.vec[i] <- computeTPR(result.table(X)$FDR, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$tpr.threshold)
        } else {
          if ('adjpvalue' %in% colnames(result.table(X))) {
            tpr.vec[i] <- computeTPR(result.table(X)$adjpvalue, 
                                     variable.annotations(X)$differential.expression,
                                     setup.parameters$tpr.threshold)
          }
        }
      }
      
      if ('nbrtpfp' %in% setup.parameters$comparisons) {
        if ('FDR' %in% colnames(result.table(X))) {
          tp.vec[i] <- computeTP(result.table(X)$FDR, 
                                 variable.annotations(X)$differential.expression,
                                 setup.parameters$nbrtpfp.threshold)
          tn.vec[i] <- computeTN(result.table(X)$FDR, 
                                 variable.annotations(X)$differential.expression,
                                 setup.parameters$nbrtpfp.threshold)
          fp.vec[i] <- computeFP(result.table(X)$FDR, 
                                 variable.annotations(X)$differential.expression,
                                 setup.parameters$nbrtpfp.threshold)
          fn.vec[i] <- computeFN(result.table(X)$FDR, 
                                 variable.annotations(X)$differential.expression,
                                 setup.parameters$nbrtpfp.threshold)
        } else {
          if ('adjpvalue' %in% colnames(result.table(X))) {
            tp.vec[i] <- computeTP(result.table(X)$adjpvalue, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$nbrtpfp.threshold)
            tn.vec[i] <- computeTN(result.table(X)$adjpvalue, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$nbrtpfp.threshold)
            fp.vec[i] <- computeFP(result.table(X)$adjpvalue, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$nbrtpfp.threshold)
            fn.vec[i] <- computeFN(result.table(X)$adjpvalue, 
                                   variable.annotations(X)$differential.expression,
                                   setup.parameters$nbrtpfp.threshold)
          }
        }
      }
      
    }
    
    if ('fracsign' %in% setup.parameters$comparisons) {
      if ('FDR' %in% colnames(result.table(X))) {
        fsn.vec[i] <- length(which(result.table(X)$FDR < 
                                     setup.parameters$fracsign.threshold))/nrow(result.table(X))
      } else {
        if ('adjpvalue' %in% colnames(result.table(X))) {
          fsn.vec[i] <- length(which(result.table(X)$adjpvalue <
                                       setup.parameters$fracsign.threshold))/nrow(result.table(X))
        }
      }
    }
    
    if ('nbrsign' %in% setup.parameters$comparisons) {
      if ('FDR' %in% colnames(result.table(X))) {
        nsn.vec[i] <- length(which(result.table(X)$FDR < 
                                     setup.parameters$fracsign.threshold))
      } else {
        if ('adjpvalue' %in% colnames(result.table(X))) {
          nsn.vec[i] <- length(which(result.table(X)$adjpvalue <
                                       setup.parameters$fracsign.threshold))
        }
      }
    }
  }
  temp.results <- NULL
  if ('typeIerror' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'typeIerror' = tIe.vec)}
  if ('fracsign' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'fracsign' = fsn.vec)}
  if ('nbrsign' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'nbrsign' = nsn.vec)}
  if ('fdr' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'fdr' = fdr.vec)}
  if ('tpr' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'tpr' = tpr.vec)}
  if ('auc' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'auc' = auc.vec)}
  if ('mcc' %in% setup.parameters$comparisons) {temp.results <- cbind(temp.results, 'mcc' = mcc.vec)}
  if ('nbrtpfp' %in% setup.parameters$comparisons) {
    temp.results <- cbind(temp.results, 'fp' = fp.vec)
    temp.results <- cbind(temp.results, 'tp' = tp.vec)
    temp.results <- cbind(temp.results, 'fn' = fn.vec)
    temp.results <- cbind(temp.results, 'tn' = tn.vec)
  }
  result.table <- data.frame(cbind(setup.parameters$file.info[, c('nbr.samples', 
                                                                  'repl', 
                                                                  'de.methods')],
                                   temp.results))
  result.table
}

# Generates the missing lines of a result table
# 
# Pads a result table generated by \code{\link{createResultTable}} by extending it with lines containing NA to replace "missing lines". This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @seealso \code{\link{createResultTable}}
# @param result.table A result table generated by \code{\link{createResultTable}}.
# @return An extended result table. This can be sent to \code{\link{plotResultTable}} to create a graphical representation in the form of a boxplot. 
# @author Charlotte Soneson
padResultTable <- function(result.table) {
  unique.spc <- unique(result.table$nbr.samples)
  unique.repl <- unique(result.table$repl)
  unique.method <- unique(result.table$de.methods)
  rt.temp.1 <- NULL
  rt.temp.2 <- NULL
  rt.temp.3 <- NULL
  for (i in seq_len(length(unique.spc))) {
    for (j in seq_len(length(unique.repl))) {
      for (k in seq_len(length(unique.method))) {
        if (length(intersect(intersect(which(result.table$nbr.samples == unique.spc[i]),
                                       which(result.table$repl == unique.repl[j])),
                             which(result.table$de.methods == unique.method[k]))) == 0) {
          rt.temp.1 <- rbind(rt.temp.1, c(unique.spc[i], unique.repl[j]))
          rt.temp.2 <- c(rt.temp.2, unique.method[k])
          rt.temp.3 <- rbind(rt.temp.3, rep(NaN, ncol(result.table) - 3))
        }
      }
    }
  }
  rt.temp <- data.frame(rt.temp.1, rt.temp.2, rt.temp.3)
  if (nrow(rt.temp) != 0) {
    colnames(rt.temp) <- colnames(result.table)
    result.table <- rbind(result.table, rt.temp)
  }
  ## Replace NAs with NaNs
  result.table[is.na(result.table)] <- NaN
  for (i in seq_len(length(unique.spc))) {
    for (j in seq_len(length(unique.method))) {
      a <- intersect(which(result.table$nbr.samples == unique.spc[i]),
                    which(result.table$de.methods == unique.method[j]))
      for (w in 4:ncol(result.table)) {
        if (all(is.na(result.table[a, w]) | result.table[a, w] == 'NaN')) {
          result.table[a, w] <- (-1000)
        }
      }
    }
  }
  result.table
}

# Plot the content of a result table in a boxplot
# 
# A function to construct a boxplot to represent a result table, obtained by \code{\link{createResultTable}} and potentially subsequently padded via \code{\link{padResultTable}}. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal).
# @param result.table R result table
# @param the.asp The measure that should be plotted. Possible values are \code{"auc"} (the area under the ROC curve), \code{"fdr"} (the false discovery rate), \code{"tpr"} (the true positive rate), \code{"typeIerror"}, \code{"mcc"} (Matthew's correlation coefficient), \code{"nbrtpfp"} (the number of TP, FP, TN, FN), \code{"nbrsign"} (the number of genes considered significant) and \code{"fracsign"} (the fraction of genes considered significant).
# @author Charlotte Soneson
plotResultTable <- function(setup.parameters, result.table, the.asp) {
  the.result <- result.table[[the.asp]]
  the.methods <- unique(result.table$de.methods)
  plot.colors <- NULL
  for (i in seq_len(length(the.methods))) {
    plot.colors <- c(plot.colors, setup.parameters$specified.colors[[the.methods[i]]])
  }
  names(plot.colors) <- the.methods
  xmax.def <- max(the.result, na.rm = TRUE)
  xmin.def <- 0
  
  if (the.asp == 'typeIerror') {
    the.xlab <- 'Type I error'
    the.threshold <- setup.parameters$typeI.threshold
    the.xmin <- setup.parameters$lower.limits$typeIerror
    the.xmax <- setup.parameters$upper.limits$typeIerror
  }
  if (the.asp == 'fdr') {
    the.xlab <- 'False discovery rate'
    the.threshold <- setup.parameters$fdr.threshold
    the.xmin <- setup.parameters$lower.limits$fdr
    the.xmax <- setup.parameters$upper.limits$fdr
  }
  if (the.asp == 'tpr') {
    the.xlab <- 'True positive rate'
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$tpr
    the.xmax <- setup.parameters$upper.limits$tpr
  }
  if (the.asp == 'mcc') {
    the.xlab <- 'MCC'
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$mcc
    the.xmax <- setup.parameters$upper.limits$mcc
  }
  if (the.asp == 'auc') {
    the.xlab <- 'AUC'
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$auc
    the.xmax <- setup.parameters$upper.limits$auc
  }
  if (the.asp == 'fracsign') {
    the.xlab <- 'Fraction significant'
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$fracsign
    the.xmax <- setup.parameters$upper.limits$fracsign
  }
  if (the.asp == "nbrsign") {
    the.xlab = "Number significant"
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$nbrsign
    the.xmax <- setup.parameters$upper.limits$nbrsign
  }
  if (the.asp == "tp") {
    the.xlab <- "Number of true positives"
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$tp
    the.xmax <- setup.parameters$upper.limits$tp
  }
  if (the.asp == "tn") {
    the.xlab = "Number of true negatives"
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$tn
    the.xmax <- setup.parameters$upper.limits$tn
  }
  if (the.asp == "fp") {
    the.xlab <- "Number of false positives"
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$fp
    the.xmax <- setup.parameters$upper.limits$fp
  }
  if (the.asp == "fn") {
    the.xlab <- "Number of false negatives"
    the.threshold <- NA
    the.xmin <- setup.parameters$lower.limits$fn
    the.xmax <- setup.parameters$upper.limits$fn
  }
  xmax <- ifelse(is.null(the.xmax), xmax.def, the.xmax)
  xmin <- ifelse(is.null(the.xmin), xmin.def, the.xmin)
  
  if (any(!is.na(the.result)) && !all(the.result == -1000)) {
    graphics::par(mar = c(5, max(nchar(the.methods)), 4, 2))
    if (any(!is.na(result.table$nbr.samples))) {
      f2.temp <- factor(paste(result.table$nbr.samples, "samples per condition"),
                        levels = paste(sort(as.numeric(unique(result.table$nbr.samples)), 
                                            decreasing = FALSE), "samples per condition"))
      result.table$nbr.samples.2 <- f2.temp
      ## The ggplot command below will throw an error: 
      ## Error in `$<-.data.frame`(`*tmp*`, "weight", value = 1) : 
      ## replacement has 1 row, data has 0
      ## if one facet is completely empty (all methods). It arises since all points of
      ## the facet are outside of the plotting range. It can be ignored.
      suppressWarnings({
        print(ggplot(result.table, aes_string(x = "de.methods", y = the.asp, 
                                              fill = "de.methods")) + 
                geom_hline(yintercept = ifelse(is.na(the.threshold), xmin - 3, the.threshold), 
                           linetype = 2) + 
                geom_boxplot(outlier.size = 0) + 
                scale_fill_manual(values = plot.colors, name = "", guide = FALSE) + 
                facet_wrap(~ nbr.samples.2, ncol = 1) + 
                geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 2) + 
                xlab("") + 
                ylab(the.xlab) + 
                theme_bw() + 
                scale_y_continuous(limits = c(xmin, xmax)) + 
                coord_flip() +
                theme(strip.text.x = element_text(size = 15), 
                      axis.text.y = element_text(size = 15),  #20
                      axis.text.x = element_text(size = 15),
                      axis.title.x = element_text(size = 15)))  #20
      })
    } else {
      suppressWarnings({
        print(ggplot(result.table, aes_string(x = "de.methods", y = the.asp, 
                                              fill = "de.methods")) +
                geom_hline(yintercept = ifelse(is.na(the.threshold), xmin - 3, the.threshold), 
                           linetype = 2) + 
                geom_boxplot(outlier.size = 0) + 
                scale_fill_manual(values = plot.colors, name = "", guide = FALSE) + 
                geom_jitter(position = position_jitter(width = 0.1, height = 0), size = 2) + 
                xlab("") + 
                ylab(the.xlab) + 
                theme_bw() + 
                scale_y_continuous(limits = c(xmin, xmax)) + 
                coord_flip() +
                theme(strip.text.x = element_text(size = 15), 
                      axis.text.y = element_text(size = 15),
                      axis.text.x = element_text(size = 15), 
                      axis.title.x = element_text(size = 15)))
      })
    }  
  } else {
    cat("No results available. No truly differentially expressed genes provided?")
  }
}

# Compute Spearman correlations among scores
# 
# Computes pairwise Spearman correlations among scores obtained from different differential expression methods applied to the same data set. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param setup.parameters List of parameters (internal). Needs to contain at least the following entries:
# \itemize{
# \item \code{incl.de.methods}
# \item \code{file.info}
# \item \code{result.table}
# }
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the correlations should be computed. The correlations are computed separately for each sample size.
# @param sel.repl The replicate number (the instance of a given simulation setting) for which the correlations should be computed.
# @author Charlotte Soneson
computeCorrelation <- function(setup.parameters, sel.nbrsamples, sel.repl) {
  all.scores <- list()
  for (i in seq_len(length(sel.nbrsamples))) {
    if (any(!is.na(sel.nbrsamples))) {
      all.scores[[sel.nbrsamples[i]]] <- list()
    } else {
      all.scores <- list()
    }
    for (j in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods",
                          setup.parameters$incl.de.methods[j])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples",
                          sel.nbrsamples[i])
      idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
      idx <- intersect(intersect(idx1, idx2), idx3)
      if (length(idx) != 0) {
        X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
        if (is.list(X)) X <- convertListTocompData(X)
        if (any(!is.na(sel.nbrsamples))) {
          all.scores[[sel.nbrsamples[i]]][[setup.parameters$incl.de.methods[j]]] <- 
            result.table(X)$score
        } else {
          all.scores[[setup.parameters$incl.de.methods[j]]] <- result.table(X)$score
        }
      }
    }
  }
  for (i in seq_len(length(sel.nbrsamples))) {
    if (any(!is.na(sel.nbrsamples))) {
      curr.methods <- names(all.scores[[sel.nbrsamples[i]]])
    } else {
      curr.methods <- names(all.scores)
    }
    Spearman.correlation <- matrix(0, length(curr.methods), length(curr.methods))
    for (j in seq_len(length(curr.methods))) {
      for (k in seq_len(length(curr.methods))) {
        if (any(!is.na(sel.nbrsamples))) {
          Spearman.correlation[j, k] <- 
            stats::cor(all.scores[[sel.nbrsamples[i]]][[curr.methods[j]]],
                       all.scores[[sel.nbrsamples[i]]][[curr.methods[k]]],
                       method = 'spearman', use = 'na.or.complete')
        } else {
          Spearman.correlation[j, k] <- 
            stats::cor(all.scores[[curr.methods[j]]],
                       all.scores[[curr.methods[k]]],
                       method = 'spearman', use = 'na.or.complete')
        }
      }
    }
    rownames(Spearman.correlation) <- colnames(Spearman.correlation) <- curr.methods
    print(Spearman.correlation)
    
    if (nrow(Spearman.correlation) > 2) {
      hc <- stats::hclust(stats::as.dist(1 - Spearman.correlation), method = 'average')
      print(levelplot(Spearman.correlation[hc$order, hc$order], xlab = '', ylab = '',
                      column.values = nrow(Spearman.correlation):1,
                      at = seq(-1, 1, length = 200),
                      scales = list(y = list(labels = curr.methods[hc$order], 
                                             at = nrow(Spearman.correlation):1),
                                    x = list(rot = 90)), 
                      col.regions = grDevices::heat.colors(200), 
                      panel = function(...) {
                        arg <- list(...)
                        panel.levelplot(...)
                        panel.text(arg$x, arg$y, round(arg$z, 2))}))
      plot(hc)
    } else {
      print(levelplot(Spearman.correlation, xlab = '', ylab = '',
                      column.values = nrow(Spearman.correlation):1,
                      at = seq(-1, 1, length = 200),
                      scales = list(y = list(labels = curr.methods, 
                                             at = nrow(Spearman.correlation):1),
                                    x = list(rot = 90)), 
                      col.regions = grDevices::heat.colors(200), 
                      panel = function(...) {
                        arg <- list(...)
                        panel.levelplot(...)
                        panel.text(arg$x, arg$y, round(arg$z, 2))}))
    }
  }
}

# Compute the overlap among sets of genes called differentially expressed by different methods
# 
# Computes the overlap among sets of genes called differentially expressed by different methods applied to the same data set. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# The regular overlap analysis (i.e., with \code{plot.type = "overlap"}) computes the overlap (the number of shared genes) between the collection of genes called differentially expressed for each pair of differential expression methods. In other words, if \eqn{A} and \eqn{B} denote the collections of differentially expressed genes obtained by two methods, the overlap is defined as \deqn{O(A,B)=|A\cap B|,}where \eqn{|\cdot|} denotes the cardinality of a set. Since the expected overlap increases with the number of genes called differentially expressed, it is also possible to perform a normalized overlap analysis, which calculates instead the Sorensen index for each pair of gene collections. The Sorensen index for the two sets \eqn{A} and \eqn{B} is defined by \deqn{S(A,B)=\frac{2|A\cap B|}{|A|+|B|}.}
# 
# @param setup.parameters List of parameters. Needs to contain at least the following entries:
# \itemize{
# \item \code{incl.de.methods}
# \item \code{file.info}
# \item \code{result.table}
# }
# @param sel.nbrsamples The sample sizes (number of samples per condition) for which the correlations should be computed. The correlations are computed separately for each sample size.
# @param sel.repl The replicate number (the instance of a given simulation setting) for which the correlations should be computed.
# @param plot.type The type of comparison that will be performed. Possible values are \code{"overlap"} and \code{"sorensen"}.
# @author Charlotte Soneson
computeOverlap <- function(setup.parameters, sel.nbrsamples, sel.repl, plot.type) {
  sign.genes <- list()
  for (i in seq_len(length(sel.nbrsamples))) {
    if (any(!is.na(sel.nbrsamples))) {
      sign.genes[[sel.nbrsamples[i]]] <- list()
    } else {
      sign.genes <- list()
    }
    for (j in seq_len(length(setup.parameters$incl.de.methods))) {
      idx1 <- findFileIdx(setup.parameters$file.info, "de.methods",
                          setup.parameters$incl.de.methods[j])
      idx2 <- findFileIdx(setup.parameters$file.info, "nbr.samples",
                          sel.nbrsamples[i])
      idx3 <- findFileIdx(setup.parameters$file.info, "repl", sel.repl)
      idx <- intersect(intersect(idx1, idx2), idx3)
      if (length(idx) != 0) {
        X <- readRDS(as.character(setup.parameters$file.info$input.files[idx]))
        if (is.list(X)) X <- convertListTocompData(X)
        if ('adjpvalue' %in% colnames(result.table(X))) {
          if (any(!is.na(sel.nbrsamples))) {
            sign.genes[[sel.nbrsamples[i]]][[setup.parameters$incl.de.methods[j]]] <-
              rownames(result.table(X))[which(result.table(X)$adjpvalue < 
                                                setup.parameters$overlap.threshold)]
          } else {
            sign.genes[[setup.parameters$incl.de.methods[j]]] <- 
              rownames(result.table(X))[which(result.table(X)$adjpvalue < 
                                                setup.parameters$overlap.threshold)]
          }
        } else {
          if ('FDR' %in% colnames(result.table(X))) {
            if (any(!is.na(sel.nbrsamples))) {
              sign.genes[[sel.nbrsamples[i]]][[setup.parameters$incl.de.methods[j]]] <- 
                rownames(result.table(X))[which(result.table(X)$FDR < 
                                                  setup.parameters$overlap.threshold)]
            } else {
              sign.genes[[setup.parameters$incl.de.methods[j]]] <- 
                rownames(result.table(X))[which(result.table(X)$FDR < 
                                                  setup.parameters$overlap.threshold)]
            }
          }
        }
      }
    }
  }
  for (i in seq_len(length(sel.nbrsamples))) {
    if (any(!is.na(sel.nbrsamples))) {
      curr.methods <- names(sign.genes[[sel.nbrsamples[i]]])
    } else {
      curr.methods <- names(sign.genes)
    }
    Sorensen.index <- matrix(0, length(curr.methods), length(curr.methods))
    Overlap <- matrix(0, length(curr.methods), length(curr.methods))
    for (j in seq_len(length(curr.methods))) {
      for (k in seq_len(length(curr.methods))) {
        if (any(!is.na(sel.nbrsamples))) {
          Sorensen.index[j, k] <- 
            computeSorensen(sign.genes[[sel.nbrsamples[i]]][[curr.methods[j]]], 
                            sign.genes[[sel.nbrsamples[i]]][[curr.methods[k]]])
          Overlap[j, k] <- length(intersect(sign.genes[[sel.nbrsamples[i]]][[curr.methods[j]]],
                                            sign.genes[[sel.nbrsamples[i]]][[curr.methods[k]]]))
        } else {
          Sorensen.index[j, k] <- computeSorensen(sign.genes[[curr.methods[j]]], 
                                                  sign.genes[[curr.methods[k]]])
          Overlap[j, k] <- length(intersect(sign.genes[[curr.methods[j]]],
                                            sign.genes[[curr.methods[k]]]))
        }
      }
    }
    rownames(Sorensen.index) <- colnames(Sorensen.index) <- curr.methods
    rownames(Overlap) <- colnames(Overlap) <- curr.methods
    if (plot.type == 'Sorensen') {
      print(Sorensen.index)
      Sorensen.noNA <- Sorensen.index
      Sorensen.noNA[is.na(Sorensen.index)] <- 0
      hc <- stats::hclust(stats::as.dist(1 - Sorensen.noNA))
      print(levelplot(Sorensen.index[hc$order, hc$order], xlab = '', ylab = '',
                      column.values = nrow(Sorensen.index):1, at = seq(0, 1, length = 100),
                      scales = list(y = list(labels=curr.methods[hc$order], 
                                             at=nrow(Sorensen.index):1),
                                    x = list(rot = 90)), 
                      col.regions = grDevices::heat.colors(100), 
                      panel = function(...) {
                        arg <- list(...)
                        panel.levelplot(...)
                        panel.text(arg$x, arg$y, round(arg$z, 2))}))
    }
    if (plot.type == 'Overlap') {
      print(Overlap)
    }
  }
}

computeSorensen <- function(x, y) {
  2 * length(intersect(x, y))/(length(x) + length(y))
}

# Shorten the names of comparison methods
# 
# A function to transform the names of comparison methods from those entered in the GUI to those used within the program. This function is generally not called by the user, the main interface for comparing differential expression methods is the \code{\link{runComparisonGUI}} function.
# 
# @param input.methods A list of comparison method identifiers, as entered in the graphical user interface.
# @return A list of transformed comparison method identifiers, as used within the package (less human readable).
# @author Charlotte Soneson
shorten.method.names <- function(input.methods) {
  transform.table <- list('AUC' = 'auc', 'TPR' = 'tpr', 'FDR' = 'fdr', 
                          'Type I error' = 'typeIerror', 'MCC' = 'mcc', 
                          'Number significant' = 'nbrsign', 
                          'Number of TP, FP, TN, FN' = 'nbrtpfp',
                          'Fraction significant' = 'fracsign', 'MA plot' = 'maplot',
                          'False discovery curves, all replicates' = 'fdcurvesall',
                          'False discovery curves, one replicate' = 'fdcurvesone',
                          'ROC, all replicates' = 'rocall', 'ROC, one replicate' = 'rocone',
                          'Overlap, one replicate' = 'overlap', 
                          'Gene score vs average expression level' = 'scorevsexpr',
                          'Sorensen index, one replicate' = 'sorensen',
                          'Spearman correlation between scores' = 'correlation', 
                          'Score distribution vs number of outliers' = 'scorevsoutlier', 
                          'FDR vs average expression level' = 'fdrvsexpr',
                          'Gene score vs signal for condition-specific genes' = 'scorevssignal')
  output.methods <- rep('', length(input.methods))
  for (i in seq_len(length(input.methods))) {
    output.methods[i] <- transform.table[[input.methods[i]]]
  }
  output.methods
}

#' Check consistency of input table to \code{\link{runComparison}}
#' 
#' Check that the \code{dataset}, \code{nbr.samples}, \code{repl} and \code{de.methods} columns of a data frame are consistent with the information provided in the input files (given in the \code{input.files} column of the data frame). If there are inconsistencies or missing information in any of the columns, replace the given information with the information in the input files. 
#' 
#' @param file.table A data frame with columns named \code{input.files} and (optionally) \code{datasets}, \code{nbr.samples}, \code{repl}, \code{de.methods}.
#' @return Returns a consistent file table defining the result files that will be used as the basis for a method comparison.
#' @author Charlotte Soneson
#' @export 
#' @examples
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 1000, 
#'                                     samples.per.cond = 5, n.diffexp = 100, 
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "voom.limma", 
#'            Rmdfunction = "voom.limma.createRmd", output.directory = tmpdir, 
#'            norm.method = "TMM")
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "edgeR.exact", 
#'            Rmdfunction = "edgeR.exact.createRmd", output.directory = tmpdir, 
#'            norm.method = "TMM",
#'            trend.method = "movingave", disp.type = "tagwise")
#' 
#' ## A correct table
#' file.table <- data.frame(input.files = file.path(tmpdir, 
#'                          c("mydata_voom.limma.rds", "mydata_edgeR.exact.rds")), 
#'                          datasets = c("mydata", "mydata"),
#'                          nbr.samples = c(5, 5),
#'                          repl = c(1, 1),
#'                          stringsAsFactors = FALSE)
#' new.table <- checkTableConsistency(file.table)
#' new.table
#' 
#' ## An incorrect table
#' file.table <- data.frame(input.files = file.path(tmpdir, 
#'                          c("mydata_voom.limma.rds", "mydata_edgeR.exact.rds")), 
#'                          datasets = c("mydata", "mydata"),
#'                          nbr.samples = c(5, 3),
#'                          repl = c(2, 1),
#'                          stringsAsFactors = FALSE)
#' new.table <- checkTableConsistency(file.table)
#' new.table
#' 
#' ## A table with missing information
#' file.table <- data.frame(input.files = file.path(tmpdir, 
#'                          c("mydata_voom.limma.rds", "mydata_edgeR.exact.rds")),
#'                          stringsAsFactors = FALSE)
#' new.table <- checkTableConsistency(file.table)
#' new.table
checkTableConsistency <- function(file.table) {
  
  if (is.null(file.table$input.files)) {
    stop("Table must have a column named input.files!")
  }
  
  datasets <- rep(NA, length(file.table$input.files))
  nbr.samples <- rep(NA, length(file.table$input.files))
  repl <- rep(NA, length(file.table$input.files))
  de.methods <- rep(NA, length(file.table$input.files))
  keep <- c()
  
  file.table <- lapply(file.table, function(x) {
    if (is.factor(x)) {
      as.character(x)
    } else {
      x
    }
  })
  
  for (i in seq_len(length(file.table$input.files))) {
    if (is.character(file.table$input.files[i]) && 
          str_detect(file.table$input.files[i], "\\.rds$")) {
      X <- readRDS(file.table$input.files[i])
      if (is.list(X)) X <- convertListTocompData(X)
      
      if (!is.logical(check_compData_results(X))) {
        message(paste("File", file.table$input.files[i], "is not a valid result object and will be removed from the table."))
      } else {
        keep <- c(keep, i)
        
        ## Check dataset (must be present in the data file)
        if (is.null(info.parameters(X)$dataset)) {
          stop(paste("'dataset' information missing in the file", file.table$input.files[i]))
        }
        if (!is.null(file.table$datasets)) {
          if (file.table$datasets[i] != info.parameters(X)$dataset) {
            message(paste("Inconsistency between provided dataset and information in file", 
                          file.table$input.files[i], ". Replacing the given information with 
                      the one given in the saved file."))
          } 
        } 
        datasets[i] <- info.parameters(X)$dataset
        
        ## Check nbr.samples
        if (is.null(info.parameters(X)$samples.per.cond)) {
          nbr.samples[i] <- NA
        } else {
          nbr.samples[i] <- info.parameters(X)$samples.per.cond
        }
        if (!is.null(file.table$nbr.samples)) {
          if (file.table$nbr.samples[i] != nbr.samples[i]) {
            message(paste("Inconsistency between provided nbr.samples and information in file", 
                          file.table$input.files[i], ". Replacing the given information with 
                      the one given in the saved file."))
          }
        }
        
        ## Check repl
        if (is.null(info.parameters(X)$repl.id)) {
          repl[i] <- NA
        } else {
          repl[i] <- info.parameters(X)$repl.id
        }
        if (!is.null(file.table$repl)) {
          if (file.table$repl[i] != repl[i]) {
            message(paste("Inconsistency between provided repl and information in file", 
                          file.table$input.files[i], ". Replacing the given information with 
                      the one given in the saved file."))
          }
        }
        
        ## Check nbr.samples
        if (is.null(method.names(X)$full.name)) {
          stop(paste("No (or misformatted) information about DE method in the file", 
                     file.table$input.files[i]))
        }
        de.methods[i] <- method.names(X)$full.name
        if (!is.null(file.table$de.methods)) {
          if (file.table$de.methods[i] != de.methods[i]) {
            message(paste("Inconsistency between provided de.method and information in file", 
                          file.table$input.files[i], ". Replacing the given information with 
                      the one given in the saved file."))
          }
        }
      }
    } else {
      message(paste(file.table$input.files[i], "is not an rds file and will be removed from the table."))
    }
  }
  ## If everything ok, add the missing info to the file.table
  file.table$datasets <- datasets
  file.table$nbr.samples <- nbr.samples
  file.table$repl <- repl
  file.table$de.methods <- de.methods
  
  file.table <- data.frame(file.table, stringsAsFactors = FALSE)
  file.table[keep, ]
}