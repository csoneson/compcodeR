#' A shiny-based GUI to the main function for running the performance
#' comparison between differential expression methods.
#'
#' This function provides a GUI to the main function for performing comparisons
#' among differential expression methods and generating a report in HTML format
#' (\code{\link{runComparison}}). It is assumed that all differential
#' expression results have been generated in advance (using e.g. the function
#' \code{\link{runDiffExp}}) and that the result \code{compData} object for
#' each data set and each differential expression method is saved separately in
#' files with the extension \code{.rds}. The function opens a graphical user
#' interface where the user can set parameter values and choose the files to be
#' used as the basis of the comparison. It is, however, possible to circumvent
#' the GUI and call the comparison function \code{\link{runComparison}} directly.
#'
#' @param input.directories A list of directories containing the result files
#'     (\code{*.rds}). All results in the provided directories will be available
#'     for inclusion in the comparison, and the selection is performed through
#'     a graphical user interface. All result objects saved in the files should
#'      be of the \code{compData} class, although list objects created by
#'      earlier versions of \code{compcodeR} are supported.
#' @param output.directory The directory where the results should be written.
#'     The subdirectory structure will be created automatically. If the
#'     directory already exists, it will be overwritten.
#' @param recursive A logical parameter indicating whether or not the search
#'     should be extended recursively to subfolders of the
#'     \code{input.directories}.
#' @param out.width The width of the figures in the final report. Will be
#'     passed on to \code{knitr} when the HTML is generated. Can be for
#'     example "800px" (see \code{knitr} documentation for more information).
#' @param upper.limits,lower.limits Lists that can be used to manually set
#'     upper and lower limits for boxplots of fdr, tpr, auc, mcc, fracsign,
#'     nbrtpfp, nbrsign and typeIerror.
#'
#' @returns
#' The function will create a comparison report, named
#' \strong{compcodeR_report<timestamp>.html}, in the \code{output.directory}.
#' It will also create subfolders named \code{compcodeR_code} and
#' \code{compcodeR_figure}, where the code used to perform the differential
#' expression analysis and the figures contained in the report, respectively,
#' will be saved. Note that if these directories already exist they will be
#' overwritten.
#'
#' @export
#' @author Charlotte Soneson
#'
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
#'
runComparisonShiny <- function(input.directories, output.directory, recursive,
                               out.width = NULL, upper.limits = NULL,
                               lower.limits = NULL) {
    ## ---------------------------------------------------------------------- ##
    ## Check input arguments
    ## ---------------------------------------------------------------------- ##
    checkClass(input.directories, "input.directories", "character")
    checkClass(output.directory, "output.directory", "character")
    checkClass(recursive, "recursive", "logical")

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

    ## List all .rds files in the input directories
    input.files.temp <- NULL
    for (input.directory in input.directories) {
        input.files.temp <- union(input.files.temp,
                                  file.path(input.directory,
                                            list.files(input.directory,
                                                       recursive = recursive,
                                                       pattern = "*.rds$")))
    }
    input.files.temp <- normalizePath(input.files.temp, winslash = "/")
    if (length(input.files.temp) == 0) {
        stop("No .rds files in the input directory")
    }

    ## ---------------------------------------------------------------------- ##
    ## Go through the files in the input.directories and see which
    ## datasets are available.
    ## ---------------------------------------------------------------------- ##
    suppressWarnings({
        all.datasets <- NULL
        all.input.files <- NULL
        for (input.file in input.files.temp) {
            temp.input <- readRDS(input.file)
            ## Check if compatible list, then convert to compData object
            if (is.list(temp.input)) {
                temp.input <- convertListTocompData(temp.input)
            }
            ## Keep only the result files (no data files or other things).
            if (is.logical(check_compData_results(temp.input))) {
                all.datasets <- c(all.datasets, info.parameters(temp.input)$dataset)
                all.input.files <- c(all.input.files, input.file)
            }
        }
        avail.datasets <- c('', unique(all.datasets))
    })

    if (length(avail.datasets) == 1) {
        stop("No result .rds files in the input directory")
    }

    ## ---------------------------------------------------------------------- ##
    ## Define shiny app
    ## ---------------------------------------------------------------------- ##
    pLayout <- shinydashboard::dashboardPage(
        skin = "purple",

        shinydashboard::dashboardHeader(
            title = paste0("compcodeR (v",
                           utils::packageVersion("compcodeR"), ")")
        ),

        shinydashboard::dashboardSidebar(
            shiny::selectInput("dataset", label = "Data set",
                               choices = avail.datasets, selected = "",
                               selectize = TRUE)
        ),

        shinydashboard::dashboardBody(
            shiny::fluidRow(
                shinydashboard::box(
                    width = 5,
                    title = "DE methods",
                    shiny::uiOutput("demethods_ui")
                ),
                shinydashboard::box(
                    width = 3,
                    title = "Replicates",
                    shiny::uiOutput("replicates_ui")
                ),
                shinydashboard::box(
                    width = 4,
                    title = "Nbr samples",
                    shiny::uiOutput("nbrsamples_ui")
                )
            ),
            shiny::fluidRow(
                shinydashboard::box(
                    width = 5,
                    title = "Parameters",
                    shiny::numericInput(
                        "adjpFDR",
                        label = "Adjusted p-value threshold for FDR",
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "adjpTPR",
                        label = "Adjusted p-value threshold for TPR",
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "adjpMCC",
                        label = "Adjusted p-value threshold for MCC",
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "nomptypeI",
                        label = "Nominal p-value threshold for type I error",
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "maxtopFDR",
                        label = paste0("Maximal number of top variables ",
                                       "in false discovery curves"),
                        value = 1500
                    ),
                    shiny::numericInput(
                        "adjpOverlap",
                        label = paste0("Adjusted p-value threshold for ",
                                       "overlap analysis"),
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "adjpFracSign",
                        label = paste0("Adjusted p-value threshold for ",
                                       "analysis of fraction significant ",
                                       "variables"),
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "adjpNbr",
                        label = paste0("Adjusted p-value threshold for ",
                                       "analysis of the number of TP, FP, ",
                                       "TN, FN variables"),
                        value = 0.05
                    ),
                    shiny::numericInput(
                        "adjpMA",
                        label = paste0("Adjusted p-value threshold for ",
                                       "coloring of variables in MA plots"),
                        value = 0.05
                    ),
                    shiny::selectInput(
                        "signalMeasure",
                        label = "Signal measure for condition-specific genes",
                        choices = c("mean", "snr"),
                        selected = "mean"
                    )
                ),
                shinydashboard::box(
                    width = 4,
                    title = "Comparisons",
                    shiny::uiOutput("comparisons_ui")
                ),
                shinydashboard::box(
                    width = 3,
                    title = "Run comparison",
                    shiny::htmlOutput("runMessage"),
                    shiny::actionButton("runcomparison", label = "Run comparison",
                                        class = "btn-primary btn-lg")
                )
            )
        )
    )

    server_function <- function(input, output, session) { # nocov start
        availableData <- shiny::reactive({
            if ((input$dataset != "")) {
                keep <- which(all.datasets == input$dataset)
                input.files <- all.input.files[keep]
                datasets <- all.datasets[keep]
                nbr.samples <- NULL
                repl <- NULL
                de.methods <- NULL
                for (input.file in input.files) {
                    temp.input <- readRDS(input.file)
                    if (is.list(temp.input)) {
                        temp.input <- convertListTocompData(temp.input)
                    }
                    ns <- info.parameters(temp.input)$samples.per.cond
                    if (is.null(ns)) {
                        ns <- NA_real_
                    }
                    nbr.samples <- c(nbr.samples, ns)

                    rp <- info.parameters(temp.input)$repl.id
                    if (is.null(rp)) {
                        rp <- NA_real_
                    }
                    repl <- c(repl, rp)

                    dm <- method.names(temp.input)$full.name
                    if (is.null(dm)) {
                        dm <- NA_character_
                    }
                    de.methods <- c(de.methods, dm)
                }
                list(
                    avail.nbr.samples = sort(as.numeric(unique(nbr.samples))),
                    avail.repl = sort(as.numeric(unique(repl))),
                    avail.de.methods = setdiff(unique(de.methods), c(NA_character_)),
                    input.files = input.files,
                    datasets = datasets,
                    nbr.samples = nbr.samples,
                    repl = repl,
                    de.methods = de.methods
                )
            } else {
                NULL
            }
        })

        output$demethods_ui <- shiny::renderUI({
            if (!is.null(availableData())) {
                shiny::checkboxGroupInput(
                    "demethods", label = "",
                    choices = availableData()$avail.de.methods,
                    selected = availableData()$avail.de.methods
                )
            } else {
                NULL
            }
        })

        output$replicates_ui <- shiny::renderUI({
            if (!is.null(availableData())) {
                shiny::checkboxGroupInput(
                    "replicates", label = "",
                    choices = availableData()$avail.repl,
                    selected = availableData()$avail.repl
                )
            } else {
                NULL
            }
        })

        output$nbrsamples_ui <- shiny::renderUI({
            if (!is.null(availableData())) {
                shiny::checkboxGroupInput(
                    "nbrsamples", label = "",
                    choices = availableData()$avail.nbr.samples,
                    selected = availableData()$avail.nbr.samples
                )
            } else {
                NULL
            }
        })

        comps <- c("AUC", "ROC, one replicate", "ROC, all replicates",
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

        output$comparisons_ui <- shiny::renderUI({
            if (!is.null(availableData())) {
                shiny::checkboxGroupInput("comparisons", label = "",
                                          choices = comps,
                                          selected = comps)
            } else {
                NULL
            }
        })

        output$runMessage <- shiny::renderText("When you are ready to run the analysis, click the button below. <br><br><strong>Important!</strong> Don't close the shiny window while the analysis is running - it will close automatically once it is done.<br><br>")

        shiny::observeEvent(input$runcomparison, {
            if(!is.null(availableData())) {
                shiny::stopApp()
                message("Be patient, your analysis is running...")

                ## Extract all parameters from the panel object
                suppressWarnings({
                    parameters <- list(incl.nbr.samples = input$nbrsamples,
                                       incl.replicates = input$replicates,
                                       incl.dataset = input$dataset,
                                       incl.de.methods = input$demethods,
                                       fdr.threshold = input$adjpFDR,
                                       tpr.threshold = input$adjpTPR,
                                       mcc.threshold = input$adjpMCC,
                                       typeI.threshold = input$nomptypeI,
                                       ma.threshold = input$adjpMA,
                                       fdc.maxvar = input$maxtopFDR,
                                       signal.measure = input$signalMeasure,
                                       overlap.threshold = input$adjpOverlap,
                                       fracsign.threshold = input$adjpFracSign,
                                       nbrtpfp.threshold = input$adjpNbr,
                                       upper.limits = upper.limits,
                                       lower.limits = lower.limits,
                                       comparisons = input$comparisons)
                    nbr.samples <- as.numeric(availableData()$nbr.samples)
                    repl <- as.numeric(availableData()$repl)
                    de.methods <- availableData()$de.methods
                    datasets <- availableData()$datasets
                    input.files <- availableData()$input.files
                })

                ## Check ranges etc for thresholds
                parameters$fdr.threshold <- checkRange(
                    parameters$fdr.threshold, "fdr.treshold", 0, 1
                )
                parameters$tpr.threshold <- checkRange(
                    parameters$tpr.threshold, "tpr.threshold", 0, 1
                )
                parameters$mcc.threshold <- checkRange(
                    parameters$mcc.threshold, "mcc.threshold", 0, 1
                )
                parameters$typeI.threshold <- checkRange(
                    parameters$typeI.threshold, "typeI.threshold", 0, 1
                )
                parameters$ma.threshold <- checkRange(
                    parameters$ma.threshold, "ma.threshold", 0, 1
                )
                parameters$overlap.threshold <- checkRange(
                    parameters$overlap.threshold, "overlap.threshold", 0, 1
                )
                parameters$fracsign.threshold <- checkRange(
                    parameters$fracsign.threshold, "fracsign.threshold", 0, 1
                )
                parameters$nbrtpfp.threshold <- checkRange(
                    parameters$nbrtpfp.threshold, "nbrtpfp.threshold", 0, 1
                )

                ## Transform the names of the comparisons to make
                parameters$comparisons <- shorten.method.names(
                    parameters$comparisons
                )

                ## Run the comparison
                file.table <- data.frame(input.files = input.files,
                                         datasets = datasets,
                                         nbr.samples = nbr.samples,
                                         repl = repl,
                                         de.methods = de.methods,
                                         stringsAsFactors = FALSE)
                shiny::withProgress({
                    runComparison(file.table, parameters,
                                  output.directory, check.table = FALSE,
                                  out.width = out.width)
                    shiny::incProgress(amount = 1)
                }, message = "Running comparison, don't close this window")
            }
        })
    } # nocov end

    shiny::shinyApp(ui = pLayout, server = server_function)
}
