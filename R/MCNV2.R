# app.R — Step 1 (inheritance via pipeline OR upload) → Step 2 (Mendelian Precision analysis)

options(shiny.maxRequestSize = 500 * 1024^2)

suppressPackageStartupMessages({
  library(shiny); library(readr); library(DT)
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
  library(patchwork)
  library(digest)  # pour clés de cache
})

tags <- shiny::tags
SIZE_LEVELS <- c("[1kb–30kb[","[10kb–30kb[","[30kb–50kb[","[50kb–100kb[",
                 "[100kb–200kb[","[200kb–500kb[","[500kb–1000kb[",">=1Mb")

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
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  raw_head <- readBin(con, what = "raw", n = n)
  is_gz <- length(raw_head) >= 2 && identical(raw_head[1:2], as.raw(c(0x1f, 0x8b)))
  if (is_gz) {
    con2 <- gzcon(file(path, "rb")); on.exit(try(close(con2), silent = TRUE), add = TRUE)
    x <- tryCatch(readChar(con2, n, useBytes = TRUE), error = function(e) "")
  } else {
    x <- tryCatch(rawToChar(raw_head), error = function(e) "")
  }
  n_tab <- lengths(regmatches(x, gregexpr("\t", x, fixed = TRUE)))
  n_com <- lengths(regmatches(x, gregexpr(",", x, fixed = TRUE)))
  if (n_tab > 0 && n_tab >= n_com) "\t" else ","
}

read_delim_smart <- function(path, delim = NULL, ...) {
  if (is.null(delim)) delim <- detect_delim(path)
  con_head <- file(path, "rb"); on.exit(close(con_head), add = TRUE)
  raw_head <- readBin(con_head, what = "raw", n = 2)
  is_gz <- length(raw_head) >= 2 && identical(raw_head[1:2], as.raw(c(0x1f, 0x8b)))
  if (is_gz) {
    con <- gzcon(file(path, "rb")); on.exit(try(close(con), silent = TRUE), add = TRUE)
    readr::read_delim(con, delim = delim, show_col_types = FALSE, progress = FALSE, ...)
  } else {
    readr::read_delim(path, delim = delim, show_col_types = FALSE, progress = FALSE, ...)
  }
}

pct <- function(x) sprintf("%0.1f%%", 100 * x)

compute_size_range <- function(size_bp) {
  cut(
    size_bp,
    breaks = c(-Inf, 10000, 30000, 50000, 100000, 200000, 500000, 1000000, Inf),
    labels = SIZE_LEVELS,
    right = FALSE
  )
}

is_nonempty_string <- function(x) is.character(x) && length(x) == 1 && nzchar(x)

apply_mapping_safely <- function(df, mapping) {
  out <- df
  for (tgt in names(mapping)) {
    src <- mapping[[tgt]]
    if (is_nonempty_string(src) && src %in% names(df)) out[[tgt]] <- df[[src]]
  }
  keep_idx <- rev(!duplicated(rev(names(out))))
  out <- out[, keep_idx, drop = FALSE]
  out
}

# ---------- UI ----------
ui <- fluidPage(
  tags$head(
    tags$style(HTML(spinner_css)),
    uiOutput("dyn_css")  # CSS dynamique pour le scroll des plots
  ),
  titlePanel("MCNV2 (Mendelian Copy Number Variation Validator)"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 uiOutput("phase_badge"),
                 conditionalPanel("output.phase == 'prep'",
                                  h4("Step 1 — Inheritance source"),
                                  radioButtons("mode", NULL,
                                               choices = c("Compute (Python pipeline)" = "calc",
                                                           "Upload an inheritance file" = "upload"),
                                               selected = "upload"
                                  ),
                                  conditionalPanel("input.mode == 'calc'",
                                                   hr(), h4("Files & parameters (Compute)"),
                                                   textInput("py_cmd", "Python", value = if (nzchar(Sys.which('python3'))) 'python3' else 'python'),
                                                   fileInput("py_script", "Python script", accept = NULL),
                                                   fileInput("cnv_tsv",  "CNV table",  accept = c(".csv",".tsv",".txt",".csv.gz",".tsv.gz",".txt.gz")),
                                                   fileInput("ped_tsv",  "Pedigree table", accept = c(".csv",".tsv",".txt",".csv.gz",".tsv.gz",".txt.gz")),
                                                   fileInput("prob_tsv", "Problematic regions", accept = c(".tsv",".tsv.gz",".txt",".txt.gz")),
                                                   radioButtons("build", "Genome build", c("38","37"), selected = "38", inline = TRUE),
                                                   numericInput("th_prob", "Problematic regions threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
                                                   numericInput("th_cnv", "Inheritance threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
                                                   hr(), h4("Execution"),
                                                   div(class = "run-status", uiOutput("run_state")),
                                                   div(class = "outdir-box", verbatimTextOutput("outdir_display")),
                                                   actionButton("run_pipe", "Run pipeline")
                                  ),
                                  conditionalPanel("input.mode == 'upload'",
                                                   hr(), h4("Upload an inheritance file"),
                                                   fileInput("inherit_file", "Inheritance file (CSV/TSV recommandé, .gz accepté)",
                                                             accept = c(".csv",".tsv",".txt",".csv.gz",".tsv.gz",".txt.gz")),
                                                   tags$small("Tip: CSV ou tab-delimited recommandé ; .gz accepté pour les gros fichiers."),
                                                   actionButton("show_inherit", "Show preview")
                                  )
                 ),
                 conditionalPanel("output.phase == 'pm'",
                                  actionButton("back_to_prep", "← Back to upload/compute", class = "btn btn-secondary"),
                                  br(), br(),
                                  h4("Step 2 — Mendelian Precision analysis"),
                                  h5("Column mapping"),
                                  uiOutput("col_map_ui"),
                                  radioButtons(
                                    "mp_mode", "Mendelian Precision mode",
                                    choices = c("CNV-level (Transmitted_CNV)" = "cnv",
                                                "Gene-level (transmitted_gene)" = "gene"),
                                    selected = "cnv", inline = TRUE
                                  ),
                                  tags$small("(CNV: dédup par CNV enfant ; Gene: exclut intergéniques)"),
                                  hr(),
                                  h4("Filters"),
                                  checkboxGroupInput("type_sel", "CNV type", choices = c("DEL","DUP"), selected = c("DEL","DUP")),
                                  uiOutput("score_ui"),
                                  uiOutput("size_ui"),
                                  uiOutput("twoalgos_ui"),
                                  uiOutput("prob_overlap_ui"),
                                  checkboxInput("exclude_qcfail", "Exclude samples that fail QC", value = TRUE),
                                  textInput("samples_exclude", "Exclude samples (IDs séparés par ,)", value = ""),
                                  hr(),
                                  h4("Gene-level filters for 4-panel"),
                                  numericInput("loeuf_cutoff", "LOEUF cutoff (≥)", value = 0.6, min = 0, max = 3, step = 0.05),
                                  fileInput("green_list", "Green genes list — supply Gene_ID (one per line)", accept = NULL),
                                  tags$small("Panels: All | LOEUF≥cutoff (gènes présents; intergéniques inclus en CNV) | No green | LOEUF≥cutoff ∧ No green."),
                                  hr(),
                                  h4("Plot options"),
                                  radioButtons("plot_type", "Type of plot:",
                                               choices = c("MP vs Score (by Size_Range)" = "mp_score",
                                                           "Global MP by Size_Range"   = "mp_size"),
                                               selected = "mp_score"
                                  ),
                                  checkboxInput("show_counts_on_points", "Show point counts on curves", value = FALSE),
                                  uiOutput("score_thresholds_ui"),
                                  tags$small("Axe MP: démarre au multiple de 5% ≤ MP min, jusqu'à 100%."),
                                  hr(),
                                  h4("Display options"),
                                  sliderInput("plot_width_px", "Plot width (px)", min = 800, max = 4000, value = 1400, step = 100),
                                  checkboxInput("plot_scroll", "Allow horizontal scroll (avoid auto-shrink)", TRUE),
                                  tags$small("Augmente la largeur et active le scroll horizontal pour voir les nombres coupés.")
                 )
    ),
    mainPanel(width = 9,
              uiOutput("main_tabs")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {
  rv <- reactiveValues(phase="prep", running=FALSE, outdir="", df=NULL, df_path=NULL)
  
  # CSS dynamique (scroll vs auto-shrink)
  output$dyn_css <- renderUI({
    if (isTRUE(input$plot_scroll)) {
      tags$style(HTML("
        #plot_del img, #plot_dup img { max-width: none; }
        .plot-scroll { overflow-x: auto; -webkit-overflow-scrolling: touch; }
      "))
    } else {
      tags$style(HTML("
        #plot_del img, #plot_dup img { max-width: 100%; }
        .plot-scroll { overflow-x: visible; }
      "))
    }
  })
  
  # Init outputs
  output$preview_tbl     <- renderDT({ datatable(data.frame()) })
  output$result_meta     <- renderText({ "" })
  output$pm_continue_btn <- renderUI({ NULL })
  output$phase           <- renderText(rv$phase)
  output$run_state       <- renderUI({ tags$span("") })
  output$outdir_display  <- renderText({ "" })
  outputOptions(output, "preview_tbl", suspendWhenHidden = FALSE)
  outputOptions(output, "result_meta", suspendWhenHidden = FALSE)
  outputOptions(output, "pm_continue_btn", suspendWhenHidden = FALSE)
  outputOptions(output, "phase", suspendWhenHidden = FALSE)
  
  output$phase_badge <- renderUI({
    tags$p(tags$b(if (rv$phase=="prep") "Step 1: Mendelian annotation" else "Step 2: Mendelian Precision analysis"))
  })
  
  output$main_tabs <- renderUI({
    if (rv$phase == "prep") {
      tabsetPanel(id = "tabs",
                  tabPanel("Mendelian annotation",
                           h4("Preview (first 50 rows)"),
                           DTOutput("preview_tbl"),
                           br(), verbatimTextOutput("result_meta"),
                           br(), uiOutput("pm_continue_btn")
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
                  tabPanel("Plot DEL",
                           div(class = "plot-scroll",
                               plotOutput("plot_del", height = 600)
                           )
                  ),
                  tabPanel("Plot DUP",
                           div(class = "plot-scroll",
                               plotOutput("plot_dup", height = 600)
                           )
                  ),
                  tabPanel("Filtered table", DTOutput("tbl")),
                  tabPanel("De novo & biomaRt",
                           fluidRow(
                             column(6, DTOutput("dnv_tbl")),
                             column(6, h5("Selection: locus & genes (biomaRt)"),
                                    verbatimTextOutput("selected_locus"),
                                    tableOutput("biomart_hits")
                             )
                           )
                  )
      )
    }
  })
  
  show_message_table <- function(txt) {
    output$preview_tbl <- renderDT({ datatable(data.frame(message=txt), options=list(dom='t')) })
    output$result_meta <- renderText({ "" })
    output$pm_continue_btn <- renderUI({ NULL })
  }
  show_message_table("Waiting for an action…")
  notify <- function(msg, duration=6) showNotification(msg, type="message", duration=duration)
  
  # -------- Upload mode --------
  observeEvent(input$show_inherit, {
    if (is.null(input$inherit_file)) { notify("No inheritance file provided."); show_message_table("No inheritance file provided."); return() }
    path <- input$inherit_file$datapath
    df <- tryCatch(read_delim_smart(path), error=function(e) NULL)
    if (is.null(df) || nrow(df) < 1) {
      notify("Unable to read the inheritance file."); show_message_table("Unable to read the inheritance file."); return()
    }
    rv$df <- df; rv$df_path <- path
    output$preview_tbl <- renderDT({
      datatable(utils::head(rv$df, 50),
                options=list(pageLength=50, lengthChange=FALSE, info=TRUE, searching=FALSE, scrollX=TRUE))
    })
    output$result_meta <- renderText({ paste0("File (upload): ", normalizePath(rv$df_path), "\nTotal rows: ", nrow(rv$df)) })
    output$pm_continue_btn <- renderUI({ actionButton("go_pm_from_btn", "Proceed to Mendelian Precision analysis", class="btn btn-primary") })
    notify("Inheritance file preview shown.")
  })
  
  # -------- Compute mode (placeholder) --------
  output$run_state <- renderUI({
    if (input$mode != "calc") return(tags$span(""))
    if (isTRUE(rv$running)) tags$span(tags$span(class="spinner"), tags$span(class="run-msg","Running pipeline…"))
    else if (nzchar(rv$outdir)) tags$span(class="run-done","Done ✅") else tags$span("")
  })
  output$outdir_display <- renderText({ if (input$mode != "calc") "" else if (nzchar(rv$outdir)) normalizePath(rv$outdir) else "" })
  observeEvent(input$run_pipe, { notify("Pipeline execution via Python not shown here to keep focus on MP logic.") })
  
  observeEvent(input$go_pm_from_btn, {
    if (is.null(rv$df) || !nrow(rv$df)) { notify("No loaded/produced table. Cannot continue."); return() }
    rv$phase <- "pm"
  })
  observeEvent(input$back_to_prep, { rv$phase <- "prep"; showNotification("Back to upload/compute.", type="message", duration=4) })
  
  # ===================== STEP 2: MP ANALYSIS =====================
  raw_cnv <- reactive({ validate(need(!is.null(rv$df) && nrow(rv$df) > 0, "No active table.")); rv$df })
  
  # Column mapping UI
  output$col_map_ui <- renderUI({
    df <- req(raw_cnv()); cols <- names(df)
    suggest <- function(cands) { cand <- intersect(cands, cols); if (length(cand)) cand[1] else "" }
    fluidRow(
      column(6, selectInput("col_geneid","Gene_ID (for green list)", choices=c("",cols), selected=suggest(c("Gene_ID","GeneID","gene_id","Gene_ID_child")))),
      column(6, selectInput("col_genename","Gene name (optional)", choices=c("",cols), selected=suggest(c("gene_name","Gene","SYMBOL","symbol")))),
      column(6, selectInput("col_type","Type", choices=c("",cols), selected=suggest(c("TYPE","type","cnv_type","TYPE_child","Type","Type_child")))),
      column(6, selectInput("col_score","Score", choices=c("",cols), selected=suggest(c("SCORE","score","MAXLOGBF","SCORE_child","SNP_child","Score")))),
      column(6, selectInput("col_sample","Sample ID", choices=c("",cols), selected=suggest(c("SampleID","sample","IID","ID","SampleID_child")))),
      column(6, selectInput("col_chr","Chromosome (CHR)", choices=c("",cols), selected=suggest(c("CHR","chrom","chromosome","CHR_child","Chr")))),
      column(6, selectInput("col_start","START (bp)", choices=c("",cols), selected=suggest(c("START","start","pos","begin","START_child","Start")))),
      column(6, selectInput("col_stop","STOP (bp)", choices=c("",cols), selected=suggest(c("STOP","stop","end","pos2","STOP_child","Stop","End")))),
      column(6, selectInput("col_loeuf","LOEUF column", choices=c("",cols), selected=suggest(c("LOEUF","loeuf","LOEUF_child")))),
      column(6, selectInput("col_trans_gene","transmitted_gene (Gene mode)", choices=c("",cols), selected=suggest(c("transmitted_gene","Transmitted_Gene")))),
      column(6, selectInput("col_trans_cnv","Transmitted_CNV (CNV mode)", choices=c("",cols), selected=suggest(c("Transmitted_CNV","transmitted_CNV","transmitted_cnv","Transmission","status","trans")))),
      column(6, selectInput("col_twoalgos","TwoAlgos column (optional)", choices=c("",cols), selected=suggest(c("TwoAlgos","TWOALGOS","two_algos","TwoAlgs","twoalgos","TwoAlgos_child"))))
    )
  })
  
  # --------- Bounds helpers ---------
  score_bounds <- reactive({
    df <- req(raw_cnv()); sc <- input$col_score
    if (!is_nonempty_string(sc) || !(sc %in% names(df))) return(NULL)
    vals <- suppressWarnings(as.numeric(df[[sc]]))
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng))) return(NULL)
    list(min = floor(rng[1]), max = ceiling(rng[2]))
  })
  size_bounds <- reactive({
    df <- req(raw_cnv()); st <- input$col_start; en <- input$col_stop
    if (!is_nonempty_string(st) || !is_nonempty_string(en) || !(st %in% names(df)) || !(en %in% names(df))) return(NULL)
    sz <- suppressWarnings(as.numeric(df[[en]]) - as.numeric(df[[st]]))
    rng <- range(sz, na.rm = TRUE)
    if (!all(is.finite(rng))) return(NULL)
    list(min = max(0, floor(rng[1])), max = max(0, ceiling(rng[2])))
  })
  
  # UI: NUMÉRIQUES pour score/size (pas de sliders)
  output$score_ui <- renderUI({
    b <- score_bounds(); if (is.null(b)) return(tags$em("Score not mapped"))
    fluidRow(
      column(6, numericInput("score_min_num", "Min score", value=b$min, min=b$min, max=b$max, step=1)),
      column(6, numericInput("score_max_num", "Max score", value=b$max, min=b$min, max=b$max, step=1))
    )
  })
  output$size_ui <- renderUI({
    b <- size_bounds(); if (is.null(b)) return(tags$em("Size not mapped (START/STOP)"))
    fluidRow(
      column(6, numericInput("size_min_num", "Min size (bp)", value=b$min, min=0, max=b$max, step=1)),
      column(6, numericInput("size_max_num", "Max size (bp)", value=b$max, min=0, max=b$max, step=1))
    )
  })
  output$twoalgos_ui <- renderUI({
    df <- req(raw_cnv()); ta <- input$col_twoalgos
    if (!is_nonempty_string(ta) || !(ta %in% names(df))) return(tags$em("TwoAlgos not mapped (optional)"))
    sliderInput("twoalgos_min","TwoAlgos ≥", min=0, max=100, value=0, step=1)
  })
  output$prob_overlap_ui <- renderUI({
    df <- req(raw_cnv())
    if (!("cnv_problematic_region_overlap" %in% names(df))) {
      return(tags$em("cnv_problematic_region_overlap not found in table"))
    }
    numericInput("prob_overlap_th", "Exclude if problematic overlap ≥", value = 0.5, min = 0, max = 1, step = 0.05)
  })
  output$score_thresholds_ui <- renderUI({
    tags$small("Seuils X: ≥min(score), puis +20 jusqu'à ≤100 (automatique).")
  })
  
  # Green list
  green_set <- reactive({
    if (is.null(input$green_list)) return(character())
    path <- input$green_list$datapath
    lines <- tryCatch(readLines(path, warn=FALSE), error=function(e) character())
    unique(trimws(lines[nzchar(trimws(lines))]))
  })
  greens_hash <- reactive({ digest(green_set()) })
  
  # Helpers
  find_end_col <- function(nms) {
    cand <- c("End","END","Stop","STOP","end","stop")
    cand[which(cand %in% nms)[1]]
  }
  add_ids <- function(d) {
    end_col <- find_end_col(names(d)); validate(need(!is.null(end_col), "Map START/STOP (end column missing)."))
    has_sample <- "sample" %in% names(d)
    d$cnv_id <- if (has_sample) paste(d$sample, d$CHR, d$START, d[[end_col]], d$type, sep="|")
    else             paste(d$CHR, d$START, d[[end_col]], d$type, sep="|")
    d$locus_id <- paste(d$CHR, d$START, d[[end_col]], d$type, sep="|")
    has_id   <- "gene_id"   %in% names(d)
    has_name <- "gene_name" %in% names(d)
    id_vec   <- if (has_id) as.character(d$gene_id) else rep(NA_character_, nrow(d))
    name_vec <- if (has_name) as.character(d$gene_name) else rep(NA_character_, nrow(d))
    d$has_gene <- (!is.na(id_vec)   & nzchar(trimws(id_vec))) |
      (!is.na(name_vec) & nzchar(trimws(name_vec)))
    d
  }
  
  # Base mapping + filtres (TwoAlgos ≥ X si mappé)
  filtered <- reactive({
    df0 <- req(raw_cnv())
    mapping <- list(
      gene_id   = input$col_geneid,
      gene_name = input$col_genename,
      type      = input$col_type,
      score     = input$col_score,
      sample    = input$col_sample,
      CHR       = input$col_chr,
      START     = input$col_start,
      STOP      = input$col_stop,
      LOEUF     = input$col_loeuf,
      transmitted_gene = input$col_trans_gene,
      transmitted_CNV  = input$col_trans_cnv,
      TwoAlgos  = input$col_twoalgos
    )
    df <- apply_mapping_safely(df0, mapping)
    
    if (all(c("START","STOP") %in% names(df))) df$size <- suppressWarnings(as.numeric(df$STOP) - as.numeric(df$START))
    if ("score" %in% names(df))    df$score    <- suppressWarnings(as.numeric(df$score))
    if ("size" %in% names(df))     df$size     <- suppressWarnings(as.numeric(df$size))
    if ("LOEUF" %in% names(df))    df$LOEUF    <- suppressWarnings(as.numeric(df$LOEUF))
    if ("TwoAlgos" %in% names(df)) df$TwoAlgos <- suppressWarnings(as.numeric(df$TwoAlgos))
    if ("size" %in% names(df))     df$Size_Range <- factor(compute_size_range(df$size), levels = SIZE_LEVELS)
    if ("type" %in% names(df))     df$type <- toupper(trimws(as.character(df$type)))
    
    if ("type" %in% names(df) && !is.null(input$type_sel)) df <- dplyr::filter(df, .data$type %in% input$type_sel)
    # Filtres bornés via numériques
    if ("score" %in% names(df) && !is.null(input$score_min_num) && !is.null(input$score_max_num)) {
      mn <- suppressWarnings(as.numeric(input$score_min_num))
      mx <- suppressWarnings(as.numeric(input$score_max_num))
      if (is.finite(mn) && is.finite(mx)) df <- dplyr::filter(df, .data$score >= min(mn,mx), .data$score <= max(mn,mx))
    }
    if ("size" %in% names(df) && !is.null(input$size_min_num) && !is.null(input$size_max_num)) {
      mn <- suppressWarnings(as.numeric(input$size_min_num))
      mx <- suppressWarnings(as.numeric(input$size_max_num))
      if (is.finite(mn) && is.finite(mx)) df <- dplyr::filter(df, .data$size >= min(mn,mx), .data$size <= max(mn,mx))
    }
    if ("TwoAlgos" %in% names(df) && !is.null(input$twoalgos_min)) {
      df <- dplyr::filter(df, !is.na(.data$TwoAlgos) & .data$TwoAlgos >= input$twoalgos_min)
    }
    # Filtre "problematic region overlap" (conserve NA, enlève >= seuil)
    if ("cnv_problematic_region_overlap" %in% names(df) && !is.null(input$prob_overlap_th)) {
      thr <- suppressWarnings(as.numeric(input$prob_overlap_th))
      if (is.finite(thr)) {
        df <- df %>% dplyr::filter(is.na(.data$cnv_problematic_region_overlap) |
                                     as.numeric(.data$cnv_problematic_region_overlap) < thr)
      }
    }
    if (!is.null(input$samples_exclude) && nzchar(input$samples_exclude) && "sample" %in% names(df)) {
      bad <- unique(trimws(unlist(strsplit(input$samples_exclude, "[,;\n\t[:space:]]+"))))
      if (length(bad)) df <- dplyr::filter(df, !(.data$sample %in% bad))
    }
    if (!("gene_name" %in% names(df)) && "gene_id" %in% names(df)) df$gene_name <- df$gene_id
    df
  }) %>% bindCache(
    rv$df_path, input$col_geneid, input$col_genename, input$col_type, input$col_score,
    input$col_sample, input$col_chr, input$col_start, input$col_stop, input$col_loeuf,
    input$col_trans_gene, input$col_trans_cnv, input$col_twoalgos,
    input$type_sel, input$score_min_num, input$score_max_num,
    input$size_min_num, input$size_max_num, input$twoalgos_min,
    input$samples_exclude, input$prob_overlap_th
  )
  
  # Version "full" avec IDs utiles
  full_with_ids <- reactive({ add_ids(filtered()) }) %>% bindCache(filtered())
  
  # ============== MP dataset ==============
  mp_df <- reactive({
    d <- filtered(); req(nrow(d) > 0)
    d <- add_ids(d)
    if (identical(input$mp_mode, "cnv")) {
      trans_col <- intersect(c("Transmitted_CNV","transmitted_CNV","transmitted_cnv","Transmission","status","trans"), names(d))[1]
      validate(need(!is.null(trans_col), "No transmission column detected for CNV mode. Map 'Transmitted_CNV' in the column mapping."))
      v <- tolower(trimws(as.character(d[[trans_col]])))
      d$trans  <- NA_real_; d$trans[v %in% c("inherited","true","t","yes","y","1")] <- 1
      d$trans[v %in% c("de novo","denovo","false","f","no","n","0")] <- 0
      d$denovo <- ifelse(is.na(d$trans), NA_real_, 1 - d$trans)
      d %>% arrange(desc(!is.na(.data$trans))) %>% distinct(.data$cnv_id, .keep_all=TRUE)
    } else {
      tg_col <- input$col_trans_gene
      validate(need(is_nonempty_string(tg_col) && tg_col %in% names(d), "Column 'transmitted_gene' required for Gene mode."))
      v <- tolower(trimws(as.character(d[[tg_col]])))
      d <- d[!(v %in% c("intergenic","", "na", "nan")) & !is.na(d[[tg_col]]), , drop=FALSE]
      keys <- intersect(c("sample","CHR","START", find_end_col(names(d)), "type","gene_name"), names(d))
      if (length(keys) >= 5) d <- d %>% distinct(across(all_of(keys)), .keep_all=TRUE)
      else if ("transcript" %in% names(d)) {
        keys2 <- intersect(c("sample","CHR","START", find_end_col(names(d)), "type","transcript"), names(d))
        d <- d %>% distinct(across(all_of(keys2)), .keep_all=TRUE)
      } else d <- d %>% distinct()
      w <- tolower(trimws(as.character(d[[tg_col]])))
      d$trans  <- NA_real_; d$trans[w %in% c("inherited","true","t","yes","y","1")] <- 1
      d$trans[w %in% c("de novo","denovo","false","f","no","n","0")] <- 0
      d$denovo <- ifelse(is.na(d$trans), NA_real_, 1 - d$trans)
      d
    }
  }) %>% bindCache(filtered(), input$mp_mode, input$col_trans_gene, input$col_trans_cnv)
  
  # --------- 4 vues (All / LOEUF / No green / Both) — inclut intergéniques en mode CNV
  views_cached <- reactive({
    greens <- green_set()
    full <- full_with_ids()
    cnv_summ <- full %>%
      dplyr::group_by(.data$cnv_id) %>%
      dplyr::summarise(
        any_gene  = any(.data$has_gene, na.rm = TRUE),
        any_green = any((.data$gene_id %in% greens | .data$gene_name %in% greens) & .data$has_gene, na.rm = TRUE),
        any_lo_lt = any(.data$has_gene & !is.na(.data$LOEUF) & (.data$LOEUF < input$loeuf_cutoff), na.rm = TRUE),
        .groups = "drop"
      )
    intergenic_ids       <- cnv_summ$cnv_id[!cnv_summ$any_gene]
    loeuf_ok_genic_ids   <- cnv_summ$cnv_id[cnv_summ$any_gene & !cnv_summ$any_lo_lt]
    green_ok_genic_ids   <- cnv_summ$cnv_id[cnv_summ$any_gene & !cnv_summ$any_green]
    
    ids_loeuf <- if (identical(input$mp_mode, "cnv")) union(intergenic_ids, loeuf_ok_genic_ids) else loeuf_ok_genic_ids
    if (length(greens)) {
      ids_green <- if (identical(input$mp_mode, "cnv")) union(intergenic_ids, green_ok_genic_ids) else green_ok_genic_ids
      ids_both  <- intersect(ids_loeuf, ids_green)
    } else {
      ids_green <- character(0)
      ids_both  <- character(0)
    }
    df_mp <- mp_df()
    list(
      all   = df_mp,
      loeuf = df_mp[df_mp$cnv_id %in% ids_loeuf, , drop = FALSE],
      green = df_mp[df_mp$cnv_id %in% ids_green, , drop = FALSE],
      both  = df_mp[df_mp$cnv_id %in% ids_both , , drop = FALSE]
    )
  }) %>% bindCache(mp_df(), greens_hash(), input$loeuf_cutoff, input$mp_mode)
  
  # ---------- Pré-agrégations ----------
  thresholds_seq <- reactive({
    req(input$plot_type == "mp_score")
    d <- mp_df(); if (!nrow(d) || !("score" %in% names(d))) return(numeric(0))
    min_sc <- suppressWarnings(as.numeric(min(d$score, na.rm = TRUE))); 
    if (!is.finite(min_sc)) return(numeric(0))
    vals <- seq(from = min_sc, to = 100, by = 20)
    if (length(vals) == 0) vals <- min_sc
    unique(vals)
  }) %>% bindCache(input$plot_type, mp_df())
  
  compute_ag_score <- function(df, ths, type_label) {
    d <- dplyr::filter(df, .data$type == type_label)
    if (!nrow(d) || !("score" %in% names(d)) || !("trans" %in% names(d))) return(tibble::tibble())
    res <- lapply(ths, function(th) {
      tmp <- d %>% dplyr::filter(.data$score >= th)
      if (!nrow(tmp)) return(NULL)
      tmp %>%
        dplyr::group_by(.data$Size_Range) %>%
        dplyr::summarise(precision = mean(.data$trans, na.rm=TRUE),
                         n = dplyr::n(), .groups="drop") %>%
        dplyr::mutate(threshold = th, TYPE = type_label)
    })
    out <- dplyr::bind_rows(res)
    if (!nrow(out)) out else dplyr::mutate(out, Size_Range = factor(.data$Size_Range, levels = SIZE_LEVELS))
  }
  
  compute_ag_size <- function(df, type_label) {
    d <- dplyr::filter(df, .data$type == type_label)
    if (!nrow(d) || !("trans" %in% names(d))) return(tibble::tibble())
    out <- d %>% dplyr::group_by(.data$Size_Range) %>%
      dplyr::summarise(total_CNVs = dplyr::n(),
                       mendelian_precision = mean(.data$trans, na.rm=TRUE),
                       .groups="drop")
    if (!nrow(out)) out else dplyr::mutate(out, Size_Range = factor(.data$Size_Range, levels = SIZE_LEVELS))
  }
  
  ag_score <- reactive({
    req(input$plot_type == "mp_score")
    ths <- thresholds_seq(); views <- views_cached()
    list(
      all   = list(DEL = compute_ag_score(views$all  , ths, "DEL"),
                   DUP = compute_ag_score(views$all  , ths, "DUP")),
      loeuf = list(DEL = compute_ag_score(views$loeuf, ths, "DEL"),
                   DUP = compute_ag_score(views$loeuf, ths, "DUP")),
      green = list(DEL = compute_ag_score(views$green, ths, "DEL"),
                   DUP = compute_ag_score(views$green, ths, "DUP")),
      both  = list(DEL = compute_ag_score(views$both , ths, "DEL"),
                   DUP = compute_ag_score(views$both , ths, "DUP"))
    )
  }) %>% bindCache(thresholds_seq(), views_cached(), input$plot_type)
  
  ag_size <- reactive({
    req(input$plot_type == "mp_size")
    views <- views_cached()
    list(
      all   = list(DEL = compute_ag_size(views$all  , "DEL"),
                   DUP = compute_ag_size(views$all  , "DUP")),
      loeuf = list(DEL = compute_ag_size(views$loeuf, "DEL"),
                   DUP = compute_ag_size(views$loeuf, "DUP")),
      green = list(DEL = compute_ag_size(views$green, "DEL"),
                   DUP = compute_ag_size(views$green, "DUP")),
      both  = list(DEL = compute_ag_size(views$both , "DEL"),
                   DUP = compute_ag_size(views$both , "DUP"))
    )
  }) %>% bindCache(views_cached(), input$plot_type)
  
  # ---------- Helpers plots ----------
  percent_labels <- function(x) paste0(round(100 * x), "%")
  mp_breaks_limits <- function(vals) {
    if (length(vals) == 0) return(list(breaks = seq(0, 1, 0.05), limits = c(0, 1)))
    m <- suppressWarnings(min(vals, na.rm = TRUE))
    if (!is.finite(m)) m <- 0
    lo <- max(0, floor((m * 100) / 5) * 5 / 100)  # arrondi vers le bas au multiple de 5%
    list(breaks = seq(lo, 1, by = 0.05), limits = c(lo, 1))
  }
  
  plot_mp_score_from_ag <- function(ag, type_label, subtitle_text, show_counts = FALSE) {
    if (is.null(ag) || !nrow(ag)) {
      return(ggplot() + theme_minimal() + ggtitle(paste("Mendelian Precision —", type_label)) +
               labs(subtitle=paste0(subtitle_text, " (no data)")))
    }
    ag$Size_Range <- factor(ag$Size_Range, levels = SIZE_LEVELS)
    br <- mp_breaks_limits(ag$precision)
    p <- ggplot(ag, aes(x = threshold, y = precision, color = Size_Range, group = Size_Range)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      scale_y_continuous(limits = br$limits, breaks = br$breaks, labels = percent_labels) +
      labs(title = paste("Mendelian Precision —", type_label),
           subtitle = subtitle_text,
           x = "Score threshold (≥)",
           y = "Mendelian Precision") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(size = 13))
    if (isTRUE(show_counts)) {
      p <- p + geom_text(aes(label = n), vjust = -0.8, size = 2.6, alpha = 0.8, check_overlap = TRUE)
    }
    p
  }
  
  plot_mp_size_from_ag <- function(ag, type_label, subtitle_text) {
    if (is.null(ag) || !nrow(ag)) {
      return(ggplot() + theme_minimal() + ggtitle(paste("Mendelian Precision —", type_label)) +
               labs(subtitle=paste0(subtitle_text, " (no data)")))
    }
    ag$Size_Range <- factor(ag$Size_Range, levels = SIZE_LEVELS)
    br <- mp_breaks_limits(ag$mendelian_precision)
    ggplot(ag, aes(x = Size_Range, y = mendelian_precision)) +
      geom_col() +
      geom_text(aes(label = paste0("n=", total_CNVs)), vjust = -0.4, size = 2.6) +
      scale_y_continuous(limits = br$limits, breaks = br$breaks, labels = percent_labels) +
      labs(title = paste("Mendelian Precision —", type_label),
           subtitle = subtitle_text,
           x = "Size_Range",
           y = "Mendelian Precision") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            plot.title  = element_text(size = 13))
  }
  
  # Sous-titres / compteurs
  counts_summary <- function(df_view, type_label) {
    if (!nrow(df_view)) return("nCNV/sample=0 | nCNV(locus)=0 | CNV×gène=0 | intergenic CNVs=0")
    endc <- find_end_col(names(df_view)); if (is.null(endc)) endc <- "STOP"
    dv <- df_view %>% dplyr::filter(.data$type == type_label)
    if (!nrow(dv)) return("nCNV/sample=0 | nCNV(locus)=0 | CNV×gène=0 | intergenic CNVs=0")
    n_by_sample <- length(unique(dv$cnv_id))
    if (!("locus_id" %in% names(dv))) dv$locus_id <- paste(dv$CHR, dv$START, dv[[endc]], dv$type, sep="|")
    n_locus <- length(unique(dv$locus_id))
    full <- full_with_ids()
    full_t <- full %>% dplyr::filter(.data$type == type_label, .data$cnv_id %in% unique(dv$cnv_id))
    n_gene_lines <- sum(full_t$has_gene, na.rm=TRUE)
    intergenic_by_cnv <- full_t %>% dplyr::group_by(.data$cnv_id) %>% dplyr::summarise(is_inter = !any(.data$has_gene), .groups="drop")
    n_inter <- sum(intergenic_by_cnv$is_inter)
    sprintf("nCNV/sample=%d | nCNV(locus)=%d | CNV×gène=%d | intergenic CNVs=%d",
            n_by_sample, n_locus, n_gene_lines, n_inter)
  }
  
  # ---- DEL plot (4 panels) ----
  output$plot_del <- renderPlot({
    views <- views_cached()
    st_all   <- counts_summary(views$all  , "DEL")
    st_loeuf <- counts_summary(views$loeuf, "DEL")
    st_green <- counts_summary(views$green, "DEL")
    st_both  <- counts_summary(views$both , "DEL")
    
    if (identical(input$plot_type, "mp_score")) {
      ag <- ag_score()
      p_all   <- plot_mp_score_from_ag(ag$all$DEL  , "DEL", paste0("All CNVs — ", st_all),   input$show_counts_on_points)
      p_lo    <- plot_mp_score_from_ag(ag$loeuf$DEL, "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf), input$show_counts_on_points)
      p_green <- plot_mp_score_from_ag(ag$green$DEL, "DEL", paste0("No green (Gene_ID) — ", st_green), input$show_counts_on_points)
      p_both  <- plot_mp_score_from_ag(ag$both$DEL , "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both), input$show_counts_on_points)
    } else {
      ag <- ag_size()
      p_all   <- plot_mp_size_from_ag(ag$all$DEL  , "DEL", paste0("All CNVs — ", st_all))
      p_lo    <- plot_mp_size_from_ag(ag$loeuf$DEL, "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf))
      p_green <- plot_mp_size_from_ag(ag$green$DEL, "DEL", paste0("No green (Gene_ID) — ", st_green))
      p_both  <- plot_mp_size_from_ag(ag$both$DEL , "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both))
    }
    (p_all | p_lo) / (p_green | p_both)
  }, width = function() input$plot_width_px)
  
  # ---- DUP plot (4 panels) ----
  output$plot_dup <- renderPlot({
    views <- views_cached()
    st_all   <- counts_summary(views$all  , "DUP")
    st_loeuf <- counts_summary(views$loeuf, "DUP")
    st_green <- counts_summary(views$green, "DUP")
    st_both  <- counts_summary(views$both , "DUP")
    
    if (identical(input$plot_type, "mp_score")) {
      ag <- ag_score()
      p_all   <- plot_mp_score_from_ag(ag$all$DUP  , "DUP", paste0("All CNVs — ", st_all),   input$show_counts_on_points)
      p_lo    <- plot_mp_score_from_ag(ag$loeuf$DUP, "DUP", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf), input$show_counts_on_points)
      p_green <- plot_mp_score_from_ag(ag$green$DUP, "DUP", paste0("No green (Gene_ID) — ", st_green), input$show_counts_on_points)
      p_both  <- plot_mp_score_from_ag(ag$both$DUP , "DUP", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both), input$show_counts_on_points)
    } else {
      ag <- ag_size()
      p_all   <- plot_mp_size_from_ag(ag$all$DUP  , "DUP", paste0("All CNVs — ", st_all))
      p_lo    <- plot_mp_size_from_ag(ag$loeuf$DUP, "DUP", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf))
      p_green <- plot_mp_size_from_ag(ag$green$DUP, "DUP", paste0("No green (Gene_ID) — ", st_green))
      p_both  <- plot_mp_size_from_ag(ag$both$DUP , "DUP", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both))
    }
    (p_all | p_lo) / (p_green | p_both)
  }, width = function() input$plot_width_px)
  
  # Table
  output$tbl <- renderDT(mp_df(), options=list(pageLength=12, scrollX=TRUE), filter="top", selection="none")
  
  # De novo table
  output$dnv_tbl <- renderDT({
    d <- mp_df()
    if (!("denovo" %in% names(d))) return(DT::datatable(tibble::tibble(message="Column 'denovo' absent"), options=list(dom='t')))
    d <- d %>% dplyr::filter(as.numeric(.data$denovo) > 0); if (!nrow(d)) return(NULL)
    datatable(d, options=list(pageLength=8, scrollX=TRUE), selection="single")
  })
  output$selected_locus <- renderText({ "" })
  output$biomart_hits <- renderTable({ NULL })
  
  # Summaries / MP (global)
  output$mp_del <- renderText({
    d <- mp_df(); if (!nrow(d) || !("type" %in% names(d))) return("NA")
    del <- dplyr::filter(d, .data$type == "DEL"); if (!nrow(del) || !("trans" %in% names(del))) return("NA")
    mp <- mean(del$trans, na.rm=TRUE); if (is.na(mp)) "NA" else pct(mp)
  })
  output$mp_dup <- renderText({
    d <- mp_df(); if (!nrow(d) || !("type" %in% names(d))) return("NA")
    dup <- dplyr::filter(d, .data$type == "DUP"); if (!nrow(dup) || !("trans" %in% names(dup))) return("NA")
    mp <- mean(dup$trans, na.rm=TRUE); if (is.na(mp)) "NA" else pct(mp)
  })
  output$n_summary <- renderText({ d <- mp_df(); paste0("Rows used for MP (mode=", input$mp_mode, "): ", nrow(d)) })
  output$trans_summary <- renderText({
    d <- mp_df(); if (!nrow(d) || !("trans" %in% names(d))) return("Transmission: not available.")
    n_inh <- sum(d$trans == 1, na.rm=TRUE); n_non <- sum(d$trans == 0, na.rm=TRUE); n_na <- sum(is.na(d$trans))
    paste0("Transmission parsed → inherited: ", n_inh, " | de novo/other: ", n_non, " | NA: ", n_na)
  })
  output$dnv_obs_exp <- renderText({ "—" })
  
  # ---------- Export ----------
  agg_for_export <- reactive({
    if (identical(input$plot_type, "mp_score")) {
      ag <- ag_score()
      make <- function(x, view) {
        bind_rows(
          dplyr::mutate(x$DEL, VIEW=view, TYPE="DEL"),
          dplyr::mutate(x$DUP, VIEW=view, TYPE="DUP")
        )
      }
      bind_rows(
        make(ag$all,   "ALL"),
        make(ag$loeuf, "LOEUF_ALLGENES_GE"),
        make(ag$green, "NO_GREEN_GENE_ID"),
        make(ag$both,  "LOEUF_ALLGENES_GE_AND_NO_GREEN")
      )
    } else {
      ag <- ag_size()
      make <- function(x, view) {
        bind_rows(
          dplyr::mutate(x$DEL, VIEW=view, TYPE="DEL"),
          dplyr::mutate(x$DUP, VIEW=view, TYPE="DUP")
        )
      }
      bind_rows(
        make(ag$all,   "ALL"),
        make(ag$loeuf, "LOEUF_ALLGENES_GE"),
        make(ag$green, "NO_GREEN_GENE_ID"),
        make(ag$both,  "LOEUF_ALLGENES_GE_AND_NO_GREEN")
      )
    }
  }) %>% bindCache(input$plot_type, views_cached(), thresholds_seq())
  
  output$dl_csv <- downloadHandler(
    filename = function() paste0("aggregated_", input$plot_type, "_", input$mp_mode, "_4panels.csv"),
    content  = function(file) { readr::write_csv(agg_for_export(), file) }
  )
  
  output$dl_png <- downloadHandler(
    filename = function() paste0("plots_", input$plot_type, "_", input$mp_mode, "_DEL_4panels.png"),
    content  = function(file) {
      views <- views_cached()
      st_all   <- counts_summary(views$all  , "DEL")
      st_loeuf <- counts_summary(views$loeuf, "DEL")
      st_green <- counts_summary(views$green, "DEL")
      st_both  <- counts_summary(views$both , "DEL")
      if (identical(input$plot_type, "mp_score")) {
        ag <- ag_score()
        p_all   <- plot_mp_score_from_ag(ag$all$DEL  , "DEL", paste0("All CNVs — ", st_all),   input$show_counts_on_points)
        p_lo    <- plot_mp_score_from_ag(ag$loeuf$DEL, "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf), input$show_counts_on_points)
        p_green <- plot_mp_score_from_ag(ag$green$DEL, "DEL", paste0("No green (Gene_ID) — ", st_green), input$show_counts_on_points)
        p_both  <- plot_mp_score_from_ag(ag$both$DEL , "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both), input$show_counts_on_points)
      } else {
        ag <- ag_size()
        p_all   <- plot_mp_size_from_ag(ag$all$DEL  , "DEL", paste0("All CNVs — ", st_all))
        p_lo    <- plot_mp_size_from_ag(ag$loeuf$DEL, "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " (all genes if any) — ", st_loeuf))
        p_green <- plot_mp_size_from_ag(ag$green$DEL, "DEL", paste0("No green (Gene_ID) — ", st_green))
        p_both  <- plot_mp_size_from_ag(ag$both$DEL , "DEL", paste0("LOEUF≥", input$loeuf_cutoff, " ∧ No green — ", st_both))
      }
      g <- (p_all | p_lo) / (p_green | p_both)
      ggsave(file, g, width=12, height=9, dpi=150)
    }
  )
}

shinyApp(ui, server)
