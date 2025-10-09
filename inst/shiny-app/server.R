#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
data("mtcars")

options(shiny.maxRequestSize = Inf)  # No upload size limit

params <- getOption("MCNV2.params", default = list())

bedtools_path <- params$bedtools_path %||% "bedtools" 
results_dir   <- params$results_dir   %||% system.file("results", package = "MCNV2")

# Define server logic required to draw a histogram
function(input, output, session) {
	
	cnv_check <- reactiveVal(FALSE)
	ped_check <- reactiveVal(FALSE)
	prb_check <- reactiveVal(FALSE)
	annot_check <- reactiveVal(FALSE)
	inherit_check <- reactiveVal(FALSE)
	preproc_check <- reactiveVal(FALSE)
	prob_regions_file <- reactiveVal(NULL)
	annot_output_file <- reactiveVal(NULL)
	inherit_output_file <- reactiveVal(NULL)
	filtered_ds <- reactiveVal(NULL) #arrow dataset to avoid loading file in memory
	ready_to_plot <- reactiveVal(FALSE)
	
	observeEvent(input$cnv_tsv, {
		req(input$cnv_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$cnv_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "cnv")
		
		cnv_check(ret$status)
		if(cnv_check()){
			output$cnv_tsv_status <- renderText(ret$msg)
		} else {
			output$cnv_tsv_status <- renderText(ret$msg)
		}
	})
	
	observeEvent(input$ped_tsv, {
		req(input$ped_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$ped_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "ped")
		
		ped_check(ret$status)
		if(ped_check()){
			output$ped_tsv_status <- renderText(ret$msg)
		} else {
			output$ped_tsv_status <- renderText(ret$msg)
		}
		
	})
	
	observeEvent(input$prb_tsv, {
		req(input$prb_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$prb_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "prb")
		
		prb_check(ret$status)
		if(prb_check()){
			output$prb_tsv_status <- renderText(ret$msg)
		} else {
			output$prb_tsv_status <- renderText(ret$msg)
		}
		
	})
	
	observeEvent(cnv_check() | ped_check(), {
		if (cnv_check() & ped_check()) {
			updateActionButton(session, "submit_preprocess", disabled = FALSE)
		} else {
			updateActionButton(session, "submit_preprocess", disabled = TRUE)
		}
	})
	
	
	observeEvent(input$submit_preprocess, {
		withProgress(message = "Running analysis...", value = 0, {
			req(input$cnv_tsv)  # Attend que le fichier soit chargé
			
			cnvs_file <- input$cnv_tsv$datapath
			
			if(prb_check()){
				prob_regions_file <- input$prb_tsv$datapath
			} else {
				prob_regions_file <- system.file("resources",
																				 paste0("problematic_regions_GRCh",input$build,".bed"), 
																				 package = "MCNV2")
			}
			
			annot_output_file(file.path(results_dir, "cnvs_annotated_by_genes.tsv"))
			
			incProgress(0.5, detail = "Annotating CNVs...")
			
			ret <- MCNV2::annotate(cnvs_file = cnvs_file, 
														 prob_regions_file = prob_regions_file, 
														 output_file = annot_output_file(), 
														 genome_version = input$build, 
														 bedtools_path = bedtools_path)
			
			if(ret == 0){
				annot_check(TRUE)
				output$annot_tsv_status <- renderText(paste("✅ Path to annotation file:\n",
																										annot_output_file()))
				output$preview_preproc_tbl <- renderDataTable({
					dat <- readr::read_tsv(annot_output_file(), n_max = 50, show_col_types = FALSE)
					datatable(dat, options = list(scrollX = TRUE), rownames = FALSE)
				})
			} else {
				annot_check(FALSE)
				output$annot_tsv_status <- renderText("❌ Gene annotation failed.\nCheck logs")
			}
			showNotification("✅ Annotation complete!", type = "message")
		})
	})
	
	observeEvent(annot_check(), {
		if (annot_check()) {
			updateActionButton(session, "submit_inheritance", disabled = FALSE)
		} else {
			updateActionButton(session, "submit_inheritance", disabled = TRUE)
		}
	})
	
	
	observeEvent(input$submit_inheritance, {
		withProgress(message = "Running analysis...", value = 0, {
			updateCollapse(session, id = "preprocess_panel", 
										 close = "Annotation table (Preview)")
			
			req(annot_output_file())  # Verifie que le fichier est valide
			req(input$ped_tsv)
			
			ped_file <- input$ped_tsv$datapath
			
			inheritance_script <- system.file("python",
																				"compute_inheritance.py", 
																				package = "MCNV2")
			
			ret <- MCNV2::compute_inheritance()
			
			incProgress(0.5, detail = "Inheritance calculation...")
			
			if(ret == 0){
				inherit_check(TRUE)
				inherit_output_file("~/workspace/projects/Support/BIOINF-219/data/DEL_DUP_inheritance_readyForShiny.tsv")
				output$inherit_tsv_status <- renderText(paste("✅ Path to inheritance file:\n",
																											inherit_output_file()))
				dat <- readr::read_tsv(inherit_output_file(), n_max = 50, show_col_types = FALSE)
				output$preview_inherit_tbl <- renderDataTable({
					datatable(dat, options = list(scrollX = TRUE), rownames = FALSE)
				})
				
			} else {
				inherit_check(FALSE)
				output$inherit_tsv_status <- renderText("❌ Inheritance calculation failed.\nCheck logs")
			}
			
			showNotification("✅ Inheritance calculation complete!", type = "message")
		})
	})
	
	observeEvent(inherit_check(), {
		if (inherit_check()) {
			updateActionButton(session, "goto_mpexploration", disabled = FALSE)
		} else {
			updateActionButton(session, "goto_mpexploration", disabled = TRUE)
		}
	})
	
	observeEvent(input$goto_mpexploration, {
		updateTabItems(session, "tabs", "mp_exploration")  # move to the "mp_exploration" tab
	})
	
	output$conditional_input <- renderUI({
		if (inherit_check()) {
			tagList(
				tags$p(paste("✅ Path to inheritance file:\n",
										 inherit_output_file()))
			)
		} else {
			tagList(
				fileInput("preproc_file", label = "Preprocessed input (mandatory)"),
				verbatimTextOutput("preproc_tsv_status"),
				selectInput("build_mp", label = "Genome build",
										choices = list("GRCh38/hg38" = 38, "GRCh37/hg19" = 37),
										selected = 38)
			)
		}
	})
	
	observeEvent(input$preproc_file, {
		req(input$preproc_file)  # Attend que le fichier soit chargé
		
		filepath <- input$preproc_file$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "preproc")
		
		preproc_check(ret$status)
		if(preproc_check()){
			inherit_output_file(filepath)
			output$preproc_tsv_status <- renderText(ret$msg)
		} else {
			output$preproc_tsv_status <- renderText(ret$msg)
		}
		
	})
	
	observeEvent(inherit_output_file(), {
		if (file.exists(inherit_output_file())) {
			
			internal_columns <- c("START_child","STOP_child", "Exon_Overlap",
														"Transcript_BP_Overlap","Size_child","LOEUF", 
														"Gnomad_Max_AF", "segmentaldup_Overlap",
														"Overlap_CNV_mother","Overlap_CNV_father",
														"overlapPR", "Size_child")
			
			spec <- readr::spec_tsv(inherit_output_file())
			
			num_cols <- names(spec$cols)[
				vapply(spec$cols, function(x) {
					inherits(x, "collector_double") || inherits(x, "collector_number")
				}, logical(1))
			]
			num_cols <- num_cols[!num_cols %in% internal_columns]
			if(length(num_cols) > 0){
				updateSelectizeInput(inputId = "quality_metric", 
														 choices = num_cols)
			} else {
				ready_to_plot(FALSE)
				# TODO: ajouter un message d'alerte!
				# TODO: desactiver le bouton submit_mpviz
			}
		} else {
			updateSelectizeInput(inputId = "quality_metric", choices = NULL)
		}
	})
	
	output$qty_metric_range_ui <- renderUI({
		req(inherit_output_file())
		req(input$quality_metric)
		
		# get range without loading the full file
		dt <- data.table::fread(inherit_output_file(), select = input$quality_metric, 
														sep = "\t", showProgress = FALSE)
		range_values <- range(dt[[input$quality_metric]], na.rm = TRUE)
		
		if(is.numeric(range_values)){
			rng_infos <- create_dynamic_ticks(range_values)
			
			sliderTextInput(
				inputId = "qty_metric_range", 
				label = paste0(input$quality_metric, " range"),
				choices = rng_infos$ticks,
				selected = c(rng_infos$min, rng_infos$max),
				grid = TRUE
			)
		} else {
			tags$p("❌ The chosen quality metric is not numeric.")
		}
		
	})
	
	observeEvent(input$submit_mpviz, {
		
		req(inherit_output_file())
		req(input$quality_metric)
		
		# convert CNV size range
		min_cnv_size <- parse_cnv_size_value(input$cnv_range[1])
		max_cnv_size <- parse_cnv_size_value(input$cnv_range[2])
		
		# get quality metric thresholds
		min_qty_metric <- as.numeric(input$qty_metric_range[1])
		max_qty_metric <- as.numeric(input$qty_metric_range[2])
		
		dataset <- arrow::open_tsv_dataset(inherit_output_file())
		dataset <- dataset %>% dplyr::filter(TYPE_child %in% input$cnv_type)
		filtered_dataset <- dataset %>% dplyr::filter(#Gnomad_Max_AF <= input$freq_cutoff,
			Exon_Overlap >= input$min_exon_ov/100,
			Transcript_BP_Overlap >= input$min_transcript_ov/100
			#Size_child >= min_cnv_size,
			#Size_child <= max_cnv_size,
			#LOEUF <= input$loeuf_cutoff
			#.data[[input$quality_metric]] >= min_qty_metric,
			#.data[[input$quality_metric]] <= max_qty_metric
			
		)

		# apply exclusion list
		if(!is.null(input$exclus_genes)){
			lines <- tryCatch(readLines(input$exclus_genes$datapath, warn = FALSE), 
												error = function(e) character())
			exclusion_list <- unique(trimws(lines[nzchar(trimws(lines))]))
			filtered_dataset <- filtered_dataset %>% dplyr::filter(!GeneID %in% exclusion_list)
		}		
		
		# update reactive filtered dataset
		filtered_ds(filtered_dataset)
		
		# calcul de la MP globale
		transmission_col <- ifelse(test = input$transmission == 1, 
															 yes = "inheritance_by_cnv", 
															 no = "inheritance_by_gene")
		mp <- filtered_ds() %>%
			dplyr::select(TYPE_child, all_of(transmission_col)) %>%
			collect() %>%
			group_by(TYPE_child, !!sym(transmission_col)) %>%
			summarise(n = n(), .groups = "drop_last") %>%
			tidyr::complete(!!sym(transmission_col) := c("de novo", "inherited"), 
											fill = list(n = 0)) %>%
			summarise(
				n_de_novo = sum(n[!!sym(transmission_col) == "de novo"]),
				n_inherited = sum(n[!!sym(transmission_col) == "inherited"]),
				prop_inherited = ifelse((n_de_novo + n_inherited) > 0, 
																n_inherited / (n_de_novo + n_inherited), 
																NA_real_),
				.groups = "drop"
			)
		
		n_cnvs_all <- dataset %>% summarize(n = n()) %>% collect()
		n_cnvs_filtered <- filtered_ds() %>% summarize(n = n()) %>% collect()
		
		output$n_filtered_cnvs <- renderValueBox({
			valueBox(
				paste0(n_cnvs_filtered$n,"/",n_cnvs_all$n), "Filtered CNVs", 
				icon = icon("filter"), color = "purple"
			)
		})
		
		output$mp_del <- renderValueBox({
			valueBox(
				ifelse(test = !is.null(mp[mp$TYPE_child == "DEL",]$prop_inherited),
							 yes = pct(round(mp[mp$TYPE_child == "DEL",]$prop_inherited, digits = 3)),
							 no = "NA"), "Global MP (DEL)", icon = icon("square-minus"), color = "red"
			)
		})
		
		output$mp_dup <- renderValueBox({
			valueBox(
				ifelse(test = !is.null(mp[mp$TYPE_child == "DUP",]$prop_inherited),
							 yes = pct(round(mp[mp$TYPE_child == "DUP",]$prop_inherited, digits = 3)),
							 no = "NA"), "Global MP (DUP)", icon = icon("square-plus"), color = "teal"
			)
		})
		
		#plot_one_view(views$all,   "DEL", input$plot_type, "All CNVs")
		output$plot_overview_del <- renderPlotly({ 
			dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(TYPE_child == "DEL"), 
																				quality_metric = input$quality_metric, 
																				transmission_col = transmission_col), 
								 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(TYPE_child == "DEL"), 
								 						 quality_metric = input$quality_metric, 
								 						 transmission_col = transmission_col))
			p <- plot_mp_vs_metric(dt = dat, 
														 title = "Mendelian Precision - DEL", 
														 subtitle = "All filtered CNVs", 
														 y_lab = "Medeliaan Precision", 
														 x_lab = paste0(input$quality_metric, " threshold (≥)"))
			
			ggplotly(p) 
		}) 
		
		output$plot_overview_dup <- renderPlotly({ 
			dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(TYPE_child == "DUP"), 
																				quality_metric = input$quality_metric, 
																				transmission_col = transmission_col), 
									 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(TYPE_child == "DUP"), 
									 						 quality_metric = input$quality_metric, 
									 						 transmission_col = transmission_col))
			p <- plot_mp_vs_metric(dt = dat, 
														 title = "Mendelian Precision - DUP", 
														 subtitle = "All filtered CNVs", 
														 y_lab = "Medeliaan Precision", 
														 x_lab = paste0(input$quality_metric, "threshold (≥)"))
			
			ggplotly(p)  
		}) 
	})
	

	
}
