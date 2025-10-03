# app.R — Step 1 (inheritance via pipeline OR upload) → Step 2 (Mendelian Precision analysis)
suppressPackageStartupMessages({
  library(shiny); library(readr); library(DT)
  library(ggplot2); library(dplyr); library(tidyr); library(stringr); library(biomaRt)
  library(patchwork)
})

tags <- shiny::tags

# ---------- Styles ----------
spinner_css <- "
.run-status { margin: 6px 0 0; min-height: 28px; }
.spinner {
  display: inline-block; width: 18px; height: 18px; border: 2px solid #cbd5e1;
  border-top-color: #0ea5e9; border-radius: 50%; animation: spin 0.8s linear infinite; vertical-align: -3px;
}
@keyframes spin { to { transform: rotate(360deg); } }
.run-msg { margin-left: 8px; color: #334155; }
.run-done { color: #2c7a7b; font-weight: 600; }
.outdir-box { padding: 8px 10px; border: 1px solid #e2e8f0; border-radius: 6px; background: #f8fafc; }
.shiny-output-error, .shiny-output-error-validation { visibility: hidden; }
"

# ---------- Utilities ----------
detect_delim <- function(path, n = 2000) {
  x <- tryCatch(readChar(path, n, useBytes = TRUE), error = function(e) "")
  n_tab <- lengths(regmatches(x, gregexpr("\t", x, fixed = TRUE)))
  n_com <- lengths(regmatches(x, gregexpr(",", x, fixed = TRUE)))
  if (n_tab > 0 && n_tab >= n_com) "\t" else ","
}
pct <- function(x) sprintf("%0.1f%%", 100 * x)

# Size_Range categories
compute_size_range <- function(size_bp) {
  cut(
    size_bp,
    breaks = c(-Inf, 10000, 30000, 50000, 100000, 200000, 500000, 1000000, Inf),
    labels = c("[1kb–10kb[", "[10kb–30kb[", "[30kb–50kb[", "[50kb–100kb[",
               "[100kb–200kb[", "[200kb–500kb[", "[500kb–1000kb[", ">=1Mb"),
    right = FALSE
  )
}

# Robust parsing of transmission columns
parse_transmission_flag <- function(x) {
  if (is.null(x)) return(list(trans = rep(NA_real_, 0), denovo = rep(NA_real_, 0)))
  v <- tolower(trimws(as.character(x)))
  num <- suppressWarnings(as.numeric(v))
  is_num <- !is.na(num)
  trans <- rep(NA_real_, length(v))
  denovo <- rep(NA_real_, length(v))
  
  # numeric 1/0
  trans[is_num]  <- ifelse(num[is_num] > 0, 1, 0)
  denovo[is_num] <- ifelse(num[is_num] > 0, 0, 1)
  
  # boolean words
  trans[!is_num & v %in% c("true","t","yes","y","inherited")]  <- 1
  trans[!is_num & v %in% c("false","f","no","n")]               <- 0
  denovo[!is_num & v %in% c("true","t","yes","y","inherited")] <- 0
  denovo[!is_num & v %in% c("false","f","no","n")]              <- 1
  
  # textual fuzzy
  trans[grepl("inher", v, perl = TRUE)] <- 1
  denovo[grepl("inher", v, perl = TRUE)] <- 0
  denovo[grepl("de\\s*_?\\s*novo|denovo", v, perl = TRUE)] <- 1
  trans[grepl("de\\s*_?\\s*novo|denovo", v, perl = TRUE)]  <- 0
  
  list(trans = trans, denovo = denovo)
}

is_nonempty_string <- function(x) is.character(x) && length(x) == 1 && nzchar(x)

apply_mapping_safely <- function(df, mapping) {
  out <- df
  for (tgt in names(mapping)) {
    src <- mapping[[tgt]]
    if (is_nonempty_string(src) && src %in% names(df)) {
      out[[tgt]] <- df[[src]]
    }
  }
  # drop duplicate names (keep last mapped)
  keep_idx <- rev(!duplicated(rev(names(out))))
  out <- out[, keep_idx, drop = FALSE]
  out
}

# ---------- UI ----------
ui <- fluidPage(
  tags$head(tags$style(HTML(spinner_css))),
  titlePanel("MCNV2 (Mendelian Copy Number Variation Validator)"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 uiOutput("phase_badge"),
                 
                 # --------- STEP 1: PREP ---------
                 conditionalPanel("output.phase == 'prep'",
                                  h4("Step 1 — Inheritance source"),
                                  radioButtons("mode", NULL,
                                               choices = c("Compute (Python pipeline)" = "calc",
                                                           "Upload an inheritance file" = "upload"),
                                               selected = "calc"),
                                  
                                  # Compute mode (no extension restrictions)
                                  conditionalPanel("input.mode == 'calc'",
                                                   hr(), h4("Files & parameters (Compute)"),
                                                   textInput("py_cmd", "Python", value = if (nzchar(Sys.which('python3'))) 'python3' else 'python'),
                                                   fileInput("py_script", "Python script (any file allowed; .py recommended)", accept = NULL),
                                                   fileInput("cnv_tsv",  "CNV table (CSV/TSV recommended)",  accept = NULL),
                                                   fileInput("ped_tsv",  "Pedigree table (CSV/TSV recommended)", accept = NULL),
                                                   fileInput("genes_tsv","Gene annotations (CSV/TSV recommended)", accept = NULL),
                                                   fileInput("prob_tsv", "Problematic regions (TSV recommended)", accept = NULL),
                                                   radioButtons("build", "Genome build", c("38","37"), selected = "38", inline = TRUE),
                                                   numericInput("th_prob", "Problematic regions threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
                                                   numericInput("th_cnv", "Inheritance threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
                                                   hr(), h4("Execution"),
                                                   div(class = "run-status", uiOutput("run_state")),
                                                   div(class = "outdir-box", verbatimTextOutput("outdir_display")),
                                                   actionButton("run_pipe", "Run pipeline")
                                  ),
                                  
                                  # Upload mode (no extension restrictions)
                                  conditionalPanel("input.mode == 'upload'",
                                                   hr(), h4("Upload an inheritance file"),
                                                   fileInput("inherit_file", "Inheritance file (CSV/TSV recommended)", accept = NULL),
                                                   tags$small("Tip: CSV or tab-delimited is recommended, but any file can be selected."),
                                                   actionButton("show_inherit", "Show preview")
                                  )
                 ),
                 
                 # --------- STEP 2: ANALYSIS ---------
                 conditionalPanel("output.phase == 'pm'",
                                  actionButton("back_to_prep", "← Back to upload/compute", class = "btn btn-secondary"),
                                  br(), br(),
                                  h4("Step 2 — Mendelian Precision analysis"),
                                  h5("Column mapping"),
                                  uiOutput("col_map_ui"),
                                  tags$small("(Tip) For biomaRt, map CHR/START/STOP if available."),
                                  hr(),
                                  
                                  h4("Filters"),
                                  checkboxGroupInput("type_sel", "CNV type", choices = c("DEL","DUP"), selected = c("DEL","DUP")),
                                  uiOutput("score_ui"),
                                  uiOutput("size_ui"),
                                  sliderInput("min_exon_ov", "% Exon overlap minimum", min = 0, max = 100, value = 0, step = 1),
                                  checkboxInput("exclude_green", "Exclude green genes (PanelApp/curation)", value = FALSE),
                                  checkboxInput("exclude_qcfail", "Exclude samples that fail QC", value = TRUE),
                                  textInput("samples_exclude", "Exclude samples (IDs separated by ,)", value = ""),
                                  
                                  hr(),
                                  h4("Gene-level filters for 4-panel"),
                                  numericInput("loeuf_cutoff", "LOEUF cutoff (≤)", value = 0.6, min = 0, max = 3, step = 0.05),
                                  fileInput("green_list", "Green genes list — supply Gene_ID (one per line)", accept = NULL),
                                  tags$small("Panels: All | LOEUF>cutoff | No Gene_ID in green list | LOEUF>cutoff ∧ No Gene_ID."),
                                  
                                  hr(),
                                  h4("Plot options"),
                                  radioButtons("plot_type", "Type of plot:",
                                               choices = c("MP vs Score (by Size_Range)" = "mp_score",
                                                           "Global MP by Size_Range"   = "mp_size"),
                                               selected = "mp_score"),
                                  uiOutput("score_thresholds_ui"),
                                  tags$small("Seuils X: ≥min(score), puis +10 jusqu'à ≥100 (automatique)."),
                                  
                                  hr(),
                                  h4("Genome build (biomaRt)"),
                                  radioButtons("genome_build", "Reference", choices = c("GRCh38","GRCh37"),
                                               selected = "GRCh38", inline = TRUE),
                                  
                                  hr(),
                                  downloadButton("dl_png", "Download figure (PNG)"),
                                  downloadButton("dl_csv", "Download aggregated data (CSV)")
                 )
    ),
    
    mainPanel(width = 9,
              uiOutput("main_tabs")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  rv <- reactiveValues(
    phase   = "prep",
    running = FALSE,
    outdir  = "",
    df      = NULL,
    df_path = NULL
  )
  
  # Phase badge
  output$phase_badge <- renderUI({
    badge <- if (rv$phase == "prep") "Step 1: Mendelian annotation" else "Step 2: Mendelian Precision analysis"
    tags$p(tags$b(badge))
  })
  output$phase <- renderText(rv$phase)
  outputOptions(output, "phase", suspendWhenHidden = FALSE)
  
  # Tabs
  output$main_tabs <- renderUI({
    if (rv$phase == "prep") {
      tabsetPanel(id = "tabs",
                  tabPanel("Mendelian annotation",
                           h4("Preview (first 50 rows)"),
                           DTOutput("preview_tbl"),
                           br(),
                           verbatimTextOutput("result_meta"),
                           br(),
                           uiOutput("pm_continue_btn")
                  )
      )
    } else {
      tabsetPanel(id = "tabs",
                  tabPanel("Overview",
                           fluidRow(
                             column(4, wellPanel(h5("Global MP (DEL)"), textOutput("mp_del"))),
                             column(4, wellPanel(h5("Global MP (DUP)"), textOutput("mp_dup"))),
                             column(4, wellPanel(h5("De novo: Observed vs Expected"), textOutput("dnv_obs_exp")))
                           ),
                           verbatimTextOutput("n_summary"),
                           verbatimTextOutput("trans_summary")
                  ),
                  tabPanel("Plot DEL", plotOutput("plot_del", height = 600)),
                  tabPanel("Plot DUP", plotOutput("plot_dup", height = 600)),
                  tabPanel("Filtered table", DTOutput("tbl")),
                  tabPanel("De novo & biomaRt",
                           fluidRow(
                             column(6, DTOutput("dnv_tbl")),
                             column(6,
                                    h5("Selection: locus & genes (biomaRt)"),
                                    verbatimTextOutput("selected_locus"),
                                    tableOutput("biomart_hits")
                             )
                           )
                  )
      )
    }
  })
  
  # Helpers (Step 1)
  show_message_table <- function(txt) {
    output$preview_tbl <- renderDT({
      datatable(data.frame(message = txt), options = list(dom = 't'))
    })
    output$result_meta <- renderText({ "" })
    output$pm_continue_btn <- renderUI({ NULL })
  }
  show_message_table("Waiting for an action…")
  notify <- function(msg, duration = 6) showNotification(msg, type = "message", duration = duration)
  
  # --------- Upload mode ---------
  observeEvent(input$show_inherit, {
    if (is.null(input$inherit_file)) {
      notify("No inheritance file provided.")
      show_message_table("No inheritance file provided.")
      return(invisible(NULL))
    }
    path <- input$inherit_file$datapath
    delim <- detect_delim(path)
    df <- tryCatch(readr::read_delim(path, delim = delim, show_col_types = FALSE, progress = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || nrow(df) < 1) {
      notify("Unable to read the inheritance file (CSV/TSV).")
      show_message_table("Unable to read the inheritance file.")
      return(invisible(NULL))
    }
    rv$df <- df; rv$df_path <- path
    
    output$preview_tbl <- renderDT({
      datatable(utils::head(rv$df, 50),
                options = list(pageLength = 50, lengthChange = FALSE, info = TRUE,
                               searching = FALSE, scrollX = TRUE))
    })
    output$result_meta <- renderText({
      paste0("File (upload): ", normalizePath(rv$df_path), "\n",
             "Total rows: ", nrow(rv$df))
    })
    output$pm_continue_btn <- renderUI({
      actionButton("go_pm_from_btn", "Proceed to Mendelian Precision analysis", class = "btn btn-primary")
    })
    notify("Inheritance file preview shown.")
  })
  
  # --------- Compute mode (Python) ---------
  output$run_state <- renderUI({
    if (input$mode != "calc") return(tags$span(""))
    if (isTRUE(rv$running)) {
      tags$span(tags$span(class = "spinner"),
                tags$span(class = "run-msg", "Running pipeline…"))
    } else if (nzchar(rv$outdir)) {
      tags$span(class = "run-done", "Done ✅")
    } else tags$span("")
  })
  output$outdir_display <- renderText({
    if (input$mode != "calc") return("")
    if (nzchar(rv$outdir)) normalizePath(rv$outdir) else ""
  })
  
  check_inputs <- function() {
    missing <- character()
    if (is.null(input$py_script)) missing <- c(missing, "Python script")
    if (is.null(input$cnv_tsv))  missing <- c(missing, "CNV table")
    if (is.null(input$ped_tsv))  missing <- c(missing, "Pedigree table")
    if (is.null(input$genes_tsv))missing <- c(missing, "Gene annotations")
    if (is.null(input$prob_tsv)) missing <- c(missing, "Problematic regions")
    if (length(missing)) {
      notify(paste0("Cannot start: missing file(s) → ", paste(missing, collapse = ", ")), duration = 8)
      return(FALSE)
    }
    TRUE
  }
  
  observeEvent(input$run_pipe, {
    if (input$mode != "calc") {
      notify("Select 'Compute (Python pipeline)' to run the pipeline.")
      return(invisible(NULL))
    }
    if (!check_inputs()) return(invisible(NULL))
    
    rv$outdir <- ""; rv$df <- NULL; rv$df_path <- NULL
    show_message_table("Running…")
    
    py <- input$py_cmd
    script <- input$py_script$datapath
    if (!file.exists(script)) {
      notify(paste0("Python script not found (upload): ", script))
      show_message_table("Python script not found.")
      return(invisible(NULL))
    }
    
    base   <- path.expand("~/Downloads")
    outdir <- file.path(base, paste0("results_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    rv$outdir <- outdir
    output$outdir_display <- renderText(normalizePath(outdir))
    
    cnv  <- input$cnv_tsv$datapath
    ped  <- input$ped_tsv$datapath
    gen  <- input$genes_tsv$datapath
    prob <- input$prob_tsv$datapath
    
    log_file <- file.path(outdir, "pipeline.stdout.log")
    err_file <- file.path(outdir, "pipeline.stderr.log")
    
    args <- c("-u", script,
              "--cnv", cnv,
              "--ped", ped,
              "--genes", gen,
              "--problematic", prob,
              "--build", input$build,
              "--problematic-threshold", as.character(input$th_prob),
              "--overlap-threshold", as.character(input$th_cnv),
              "--outdir", outdir)
    
    rv$running <- TRUE
    withProgress(message = "Running pipeline…", value = 0, {
      incProgress(0.2, detail = "Init")
      status <- tryCatch(
        system2(py, args, stdout = log_file, stderr = err_file, wait = TRUE),
        error = function(e) 1L
      )
      incProgress(0.6, detail = "Reading final TSV")
      
      final_path <- file.path(outdir, "children_with_transmission.tsv")
      if (!file.exists(final_path)) {
        cand <- list.files(outdir, pattern = "(?i)\\.tsv$", full.names = TRUE, recursive = TRUE)
        if (length(cand)) {
          cand <- cand[order(file.info(cand)$mtime, decreasing = TRUE)]
          final_path <- cand[1]
        }
      }
      
      if (!file.exists(final_path)) {
        incProgress(0.2, detail = "No TSV detected")
        rv$running <- FALSE
        notify(paste0("No readable final TSV found in: ", normalizePath(outdir)), duration = 10)
        show_message_table("No final TSV detected in the output directory.")
        return(invisible(NULL))
      }
      
      df <- tryCatch(readr::read_tsv(final_path, show_col_types = FALSE, progress = FALSE),
                     error = function(e) NULL)
      if (is.null(df) || nrow(df) < 1) {
        incProgress(0.2, detail = "Failed to read TSV")
        rv$running <- FALSE
        notify(paste0("Unable to read TSV: ", normalizePath(final_path)), duration = 10)
        show_message_table("Unable to read the final TSV.")
        return(invisible(NULL))
      }
      
      rv$df <- df
      rv$df_path <- final_path
      
      output$preview_tbl <- renderDT({
        datatable(utils::head(rv$df, 50),
                  options = list(pageLength = 50, lengthChange = FALSE, info = TRUE,
                                 searching = FALSE, scrollX = TRUE))
      })
      output$result_meta <- renderText({
        paste0("File (compute): ", normalizePath(rv$df_path), "\n",
               "Total rows: ", nrow(rv$df))
      })
      output$pm_continue_btn <- renderUI({
        actionButton("go_pm_from_btn", "Proceed to Mendelian Precision analysis", class = "btn btn-primary")
      })
      
      incProgress(0.2, detail = "Done")
    })
    rv$running <- FALSE
    
    showNotification(paste("Pipeline finished. Results:", normalizePath(rv$outdir)),
                     type = "message", duration = 8)
  })
  
  # Continue to PM analysis
  observeEvent(input$go_pm_from_btn, {
    if (is.null(rv$df) || !nrow(rv$df)) {
      notify("No loaded/produced table. Cannot continue.")
      return(invisible(NULL))
    }
    rv$phase <- "pm"
  })
  
  # Back button
  observeEvent(input$back_to_prep, {
    rv$phase <- "prep"
    showNotification("Back to upload/compute.", type = "message", duration = 4)
  })
  
  # ===================== STEP 2: MP ANALYSIS =====================
  raw_cnv <- reactive({
    validate(need(!is.null(rv$df) && nrow(rv$df) > 0, "No active table."))
    rv$df
  })
  
  # Column mapping UI
  output$col_map_ui <- renderUI({
    df <- req(raw_cnv()); cols <- names(df)
    suggest <- function(cands) { cand <- intersect(cands, cols); if (length(cand)) cand[1] else "" }
    fluidRow(
      column(6, selectInput("col_geneid", "Gene_ID (for green list)", choices = c("", cols),
                            selected = suggest(c("Gene_ID","GeneID","gene_id","Gene_ID_child")))),
      column(6, selectInput("col_type",  "Type", choices = c("", cols),
                            selected = suggest(c("TYPE","type","cnv_type","TYPE_child","Type_child")))),
      column(6, selectInput("col_score", "Score", choices = c("", cols),
                            selected = suggest(c("SCORE","score","MAXLOGBF","SCORE_child","SNP_child")))),
      column(6, selectInput("col_sample","Sample ID", choices = c("", cols),
                            selected = suggest(c("SampleID","sample","IID","ID","SampleID_child")))),
      column(6, selectInput("col_chr",   "Chromosome (CHR)", choices = c("", cols),
                            selected = suggest(c("CHR","chrom","chromosome","CHR_child")))),
      column(6, selectInput("col_start", "START (bp)", choices = c("", cols),
                            selected = suggest(c("START","start","pos","begin","START_child")))),
      column(6, selectInput("col_stop",  "STOP (bp)", choices = c("", cols),
                            selected = suggest(c("STOP","stop","end","pos2","STOP_child")))),
      column(6, selectInput("col_transcnv","Transmission", choices = c("", cols),
                            selected = suggest(c("transmitted_CNV","is_transmitted","status","transmitted_gene","is_gene_transmitted")))),
      column(6, selectInput("col_loeuf","LOEUF column", choices = c("", cols),
                            selected = suggest(c("LOEUF","loeuf","LOEUF_child")))
      )
    )
  })
  
  # Score slider, size slider, thresholds "UI" text
  output$score_ui <- renderUI({
    df <- req(raw_cnv()); sc <- input$col_score
    if (!is_nonempty_string(sc) || !(sc %in% names(df))) return(tags$em("Score not mapped"))
    vals <- suppressWarnings(as.numeric(df[[sc]]))
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng))) return(tags$em("Score is not numeric"))
    sliderInput("score_rng", "Filter score", min = floor(rng[1]), max = ceiling(rng[2]), value = rng,
                step = max(1, (rng[2]-rng[1])/100))
  })
  output$size_ui <- renderUI({
    df <- req(raw_cnv()); st <- input$col_start; en <- input$col_stop
    if (!is_nonempty_string(st) || !is_nonempty_string(en) || !(st %in% names(df)) || !(en %in% names(df))) {
      return(tags$em("Size not mapped (START/STOP)"))
    }
    sz <- suppressWarnings(as.numeric(df[[en]]) - as.numeric(df[[st]]))
    rng <- range(sz, na.rm = TRUE)
    sliderInput("size_rng", "Filter size (bp)", min = 0, max = max(ceiling(rng[2]), 1), value = rng,
                step = max(1, round((rng[2]-rng[1])/100)))
  })
  output$score_thresholds_ui <- renderUI({
    tags$small("Seuils utilisés: ≥min(score), puis +10 jusqu'à ≥100 (automatique).")
  })
  
  # Read green list (Gene_ID)
  green_set <- reactive({
    if (is.null(input$green_list)) return(character())
    path <- input$green_list$datapath
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    unique(trimws(lines[nzchar(trimws(lines))]))
  })
  
  # Base mapping + derived columns and filters
  filtered <- reactive({
    df0 <- req(raw_cnv())
    # mapping
    mapping <- list(
      gene_id = input$col_geneid,
      type    = input$col_type,
      score   = input$col_score,
      sample  = input$col_sample,
      CHR     = input$col_chr,
      START   = input$col_start,
      STOP    = input$col_stop
    )
    # choose a transmission source
    trans_src <- NULL
    for (cand in c(input$col_transcnv, "transmitted_CNV", "transmitted_gene", "status", "is_transmitted")) {
      if (is_nonempty_string(cand) && cand %in% names(df0)) { trans_src <- cand; break }
    }
    if (!is.null(trans_src)) mapping$transmitted_CNV <- trans_src
    
    df <- apply_mapping_safely(df0, mapping)
    
    # numeric fields
    if (all(c("START","STOP") %in% names(df))) {
      df$size <- suppressWarnings(as.numeric(df$STOP) - as.numeric(df$START))
    }
    if ("score" %in% names(df)) df$score <- suppressWarnings(as.numeric(df$score))
    if ("size"  %in% names(df)) df$size  <- suppressWarnings(as.numeric(df$size))
    
    # LOEUF column mapped as 'LOEUF'
    if (is_nonempty_string(input$col_loeuf) && input$col_loeuf %in% names(df0)) {
      df$LOEUF <- suppressWarnings(as.numeric(df0[[input$col_loeuf]]))
    }
    
    # transmission flags
    if ("transmitted_CNV" %in% names(df)) {
      pv <- parse_transmission_flag(df$transmitted_CNV)
      df$trans  <- pv$trans
      df$denovo <- pv$denovo
    }
    
    # Size_Range
    if ("size" %in% names(df)) {
      df$Size_Range <- compute_size_range(df$size)
    }
    
    # global filters
    if ("type"  %in% names(df) && !is.null(input$type_sel)) df <- dplyr::filter(df, .data$type %in% input$type_sel)
    if ("score" %in% names(df) && !is.null(input$score_rng)) df <- dplyr::filter(df, .data$score >= input$score_rng[1], .data$score <= input$score_rng[2])
    if ("size"  %in% names(df) && !is.null(input$size_rng))  df <- dplyr::filter(df, .data$size  >= input$size_rng[1],  .data$size  <= input$size_rng[2])
    
    if (!is.null(input$samples_exclude) && nzchar(input$samples_exclude) && "sample" %in% names(df)) {
      bad <- unique(trimws(unlist(strsplit(input$samples_exclude, "[,;\n\t[:space:]]+"))))
      if (length(bad)) df <- dplyr::filter(df, !(.data$sample %in% bad))
    }
    df
  })
  
  # Summaries / MP (global)
  output$mp_del <- renderText({
    df <- filtered(); if (!nrow(df) || !("type" %in% names(df))) return("NA")
    del <- dplyr::filter(df, .data$type == "DEL"); if (!nrow(del) || !("trans" %in% names(del))) return("NA")
    mp <- mean(del$trans, na.rm = TRUE)
    if (is.na(mp)) "NA" else pct(mp)
  })
  output$mp_dup <- renderText({
    df <- filtered(); if (!nrow(df) || !("type" %in% names(df))) return("NA")
    dup <- dplyr::filter(df, .data$type == "DUP"); if (!nrow(dup) || !("trans" %in% names(dup))) return("NA")
    mp <- mean(dup$trans, na.rm = TRUE)
    if (is.na(mp)) "NA" else pct(mp)
  })
  output$n_summary <- renderText({
    df <- filtered(); paste0("Filtered CNVs: ", nrow(df))
  })
  output$trans_summary <- renderText({
    df <- filtered()
    if (!nrow(df) || !("trans" %in% names(df))) return("Transmission: not available.")
    n_inh <- sum(df$trans == 1, na.rm = TRUE)
    n_non <- sum(df$trans == 0, na.rm = TRUE)
    n_na  <- sum(is.na(df$trans))
    paste0("Transmission parsed → inherited: ", n_inh, " | de novo/other: ", n_non, " | NA: ", n_na)
  })
  output$dnv_obs_exp <- renderText({ "—" })
  
  # ===== Aggregations for plots =====
  # Seuils fixes pour l'axe X: >=minScore, puis +10 jusqu'à >=100
  thresholds_seq <- reactive({
    req(input$plot_type == "mp_score")
    min_sc <- suppressWarnings(as.numeric(req(input$score_rng)[1]))
    if (!is.finite(min_sc)) return(numeric(0))
    ths <- c(min_sc, seq(from = min_sc + 10, to = 100, by = 10))
    unique(ths)
  })
  
  mp_vs_score_by_size <- function(df, type_label, thresholds_vec) {
    d <- dplyr::filter(df, .data$type == type_label)
    if (!nrow(d) || !("score" %in% names(d)) || !("trans" %in% names(d))) return(tibble::tibble())
    res <- lapply(thresholds_vec, function(th) {
      tmp <- d %>% dplyr::filter(.data$score >= th)
      if (!nrow(tmp)) return(NULL)
      tmp %>%
        dplyr::group_by(.data$Size_Range) %>%
        dplyr::summarise(precision = mean(.data$trans, na.rm = TRUE),
                         n = dplyr::n(),
                         .groups = "drop") %>%
        dplyr::mutate(threshold = th, TYPE = type_label)
    })
    dplyr::bind_rows(res)
  }
  
  mp_by_size <- function(df, type_label) {
    d <- dplyr::filter(df, .data$type == type_label)
    if (!nrow(d) || !("trans" %in% names(d))) return(tibble::tibble())
    d %>%
      dplyr::group_by(.data$Size_Range) %>%
      dplyr::summarise(
        total_CNVs = dplyr::n(),
        mendelian_precision = mean(.data$trans, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # ---- Build four filtered views (All | LOEUF>cutoff | No green (by Gene_ID) | LOEUF>cutoff ∧ No green) ----
  make_four_views <- function(df) {
    greens <- green_set()
    have_geneid <- "gene_id" %in% names(df)
    have_loeuf  <- "LOEUF" %in% names(df)
    cutv <- input$loeuf_cutoff
    
    list(
      all   = df,
      loeuf = if (have_loeuf) dplyr::filter(df, !is.na(.data$LOEUF) & .data$LOEUF > cutv) else df[0,],
      green = if (have_geneid && length(greens))
        dplyr::filter(df, !(.data$gene_id %in% greens)) else df[0,],
      both  = if (have_geneid && have_loeuf && length(greens))
        dplyr::filter(df, !is.na(.data$LOEUF) & .data$LOEUF > cutv & !(.data$gene_id %in% greens))
      else df[0,]
    )
  }
  
  # ---- One panel (single subplot) for a given view ----
  plot_one_view <- function(df_view, type_label, mode, subtitle_tag) {
    if (!nrow(df_view) || !all(c("type","trans","Size_Range") %in% names(df_view))) {
      return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
    }
    if (identical(mode, "mp_score")) {
      ths <- thresholds_seq()
      ag <- mp_vs_score_by_size(df_view, type_label, ths)
      if (!nrow(ag)) return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
      ag$Size_Range <- factor(ag$Size_Range,
                              levels = c("[1kb–10kb[","[10kb–30kb[","[30kb–50kb[","[50kb–100kb[",
                                         "[100kb–200kb[","[200kb–500kb[","[500kb–1000kb[",">=1Mb"))
      ggplot(ag, aes(x = threshold, y = precision, color = Size_Range)) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.8) +
        geom_text(aes(label = n), vjust = -0.8, size = 2.6, alpha = 0.8, check_overlap = TRUE) +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
        labs(title = paste("Mendelian Precision —", type_label),
             subtitle = subtitle_tag,
             x = "Score threshold (≥)", y = "Mendelian Precision") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(size = 13))
    } else {
      ag <- mp_by_size(df_view, type_label)
      if (!nrow(ag)) return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
      ag$Size_Range <- factor(ag$Size_Range,
                              levels = c("[1kb–10kb[","[10kb–30kb[","[30kb–50kb[","[50kb–100kb[",
                                         "[100kb–200kb[","[200kb–500kb[","[500kb–1000kb[",">=1Mb"))
      ggplot(ag, aes(x = Size_Range, y = mendelian_precision)) +
        geom_col() +
        geom_text(aes(label = paste0("n=", total_CNVs)), vjust = -0.4, size = 2.6) +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
        labs(title = paste("Mendelian Precision —", type_label),
             subtitle = subtitle_tag,
             x = "Size_Range", y = "Mendelian Precision") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 25, hjust = 1),
              plot.title = element_text(size = 13))
    }
  }
  
  # ===== DEL plot (4 panels) =====
  output$plot_del <- renderPlot({
    df <- filtered()
    validate(need(nrow(df) > 0, "No data"))
    views <- make_four_views(df)
    
    p_all   <- plot_one_view(views$all,   "DEL", input$plot_type, "All CNVs")
    p_lo    <- plot_one_view(views$loeuf, "DEL", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff))
    p_green <- plot_one_view(views$green, "DEL", input$plot_type, "No green (by Gene_ID)")
    p_both  <- plot_one_view(views$both,  "DEL", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff, " ∧ No green (Gene_ID)"))
    
    (p_all | p_lo) / (p_green | p_both)
  })
  
  # ===== DUP plot (4 panels) =====
  output$plot_dup <- renderPlot({
    df <- filtered()
    validate(need(nrow(df) > 0, "No data"))
    views <- make_four_views(df)
    
    p_all   <- plot_one_view(views$all,   "DUP", input$plot_type, "All CNVs")
    p_lo    <- plot_one_view(views$loeuf, "DUP", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff))
    p_green <- plot_one_view(views$green, "DUP", input$plot_type, "No green (by Gene_ID)")
    p_both  <- plot_one_view(views$both,  "DUP", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff, " ∧ No green (Gene_ID)"))
    
    (p_all | p_lo) / (p_green | p_both)
  })
  
  # Table
  output$tbl <- renderDT(filtered(), options = list(pageLength = 12, scrollX = TRUE), filter = "top", selection = "none")
  
  # De novo & biomaRt (shown only if 'denovo' exists)
  output$dnv_tbl <- renderDT({
    df <- filtered()
    if (!("denovo" %in% names(df))) return(DT::datatable(tibble::tibble(message="Column 'denovo' absent"), options=list(dom='t')))
    d <- df %>% dplyr::filter(as.numeric(.data$denovo) > 0); if (!nrow(d)) return(NULL)
    datatable(d, options = list(pageLength = 8, scrollX = TRUE), selection = "single")
  })
  output$selected_locus <- renderText({ "" })
  output$biomart_hits <- renderTable({ NULL })
  
  # Downloads (aggregated data / plots)
  agg_for_export <- reactive({
    df <- filtered()
    if (!nrow(df)) return(tibble::tibble())
    if (identical(input$plot_type, "mp_score")) {
      ths <- req(thresholds_seq())
      # export both types, all four views concatenated with tags
      views <- make_four_views(df)
      lab <- function(k) c(all="ALL", loeuf="LOEUF_GT", green="NO_GREEN_GENE_ID", both="LOEUF_GT_AND_NO_GREEN_GENE_ID")[[k]]
      bind <- function(v, typ) if (!nrow(v)) tibble::tibble() else {
        dplyr::bind_rows(
          mp_vs_score_by_size(v, "DEL", ths) %>% dplyr::mutate(TYPE="DEL", VIEW=typ),
          mp_vs_score_by_size(v, "DUP", ths) %>% dplyr::mutate(TYPE="DUP", VIEW=typ)
        )
      }
      dplyr::bind_rows(
        bind(views$all,   lab("all")),
        bind(views$loeuf, lab("loeuf")),
        bind(views$green, lab("green")),
        bind(views$both,  lab("both"))
      )
    } else {
      views <- make_four_views(df)
      lab <- function(k) c(all="ALL", loeuf="LOEUF_GT", green="NO_GREEN_GENE_ID", both="LOEUF_GT_AND_NO_GREEN_GENE_ID")[[k]]
      bind <- function(v, typ) if (!nrow(v)) tibble::tibble() else {
        dplyr::bind_rows(
          mp_by_size(v, "DEL") %>% dplyr::mutate(TYPE="DEL", VIEW=typ),
          mp_by_size(v, "DUP") %>% dplyr::mutate(TYPE="DUP", VIEW=typ)
        )
      }
      dplyr::bind_rows(
        bind(views$all,   lab("all")),
        bind(views$loeuf, lab("loeuf")),
        bind(views$green, lab("green")),
        bind(views$both,  lab("both"))
      )
    }
  })
  
  output$dl_csv <- downloadHandler(
    filename = function() paste0("aggregated_", input$plot_type, "_4panels.csv"),
    content = function(file) { readr::write_csv(agg_for_export(), file) }
  )
  output$dl_png <- downloadHandler(
    filename = function() paste0("plots_", input$plot_type, "_DEL_4panels.png"),
    content = function(file) {
      df <- filtered()
      views <- make_four_views(df)
      p_all   <- plot_one_view(views$all,   "DEL", input$plot_type, "All CNVs")
      p_lo    <- plot_one_view(views$loeuf, "DEL", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff))
      p_green <- plot_one_view(views$green, "DEL", input$plot_type, "No green (by Gene_ID)")
      p_both  <- plot_one_view(views$both,  "DEL", input$plot_type, paste0("LOEUF > ", input$loeuf_cutoff, " ∧ No green (Gene_ID)"))
      g <- (p_all | p_lo) / (p_green | p_both)
      ggsave(file, g, width = 12, height = 9, dpi = 150)
    }
  )
}

shinyApp(ui, server)
