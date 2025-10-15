#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(datamods)

options(shiny.maxRequestSize = Inf)  # No upload size limit
options(arrow.pull_as_vector = TRUE)

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
	genes_check <- reactiveVal(FALSE)
	
	prob_regions_file <- reactiveVal(NULL)
	annot_output_file <- reactiveVal(NULL)
	inherit_output_file <- reactiveVal(NULL)
	exclusion_list <- reactiveVal(NULL)
	
	quality_metrics <- reactiveVal(NULL)
	filtered_ds <- reactiveVal(NULL) #arrow dataset to avoid loading file in memory
	
	baseline_DEL_plot <-  reactiveVal(NULL)
	baseline_DUP_plot <-  reactiveVal(NULL)
	
	observeEvent(input$cnv_tsv, {
		req(input$cnv_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$cnv_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "cnv")
		cnv_check(ret$status)
		output$cnv_tsv_status <- renderText(ret$msg)
	})
	
	observeEvent(input$ped_tsv, {
		req(input$ped_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$ped_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "ped")
		ped_check(ret$status)
		output$ped_tsv_status <- renderText(ret$msg)
		
	})
	
	observeEvent(input$prb_tsv, {
		req(input$prb_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$prb_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, file_type = "prb")
		prb_check(ret$status)
		output$prb_tsv_status <- renderText(ret$msg)
		
	})
	
	observeEvent(input$exclus_genes, {
		req(input$exclus_genes)  # Attend que le fichier soit chargé
		
		filepath <- input$exclus_genes$datapath
		
		ret <- load_gene_ids(filepath)
		genes_check(ret$status)
		output$exclus_genes_status <- renderText(ret$msg)
		if(genes_check()){
			exclusion_list(ret$genes)
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
																				 paste0("problematic_regions_GRCh",
																				 			 input$build,".bed"), 
																				 package = "MCNV2")
			}
			
			annot_output_file(file.path(results_dir, paste0("cnvs_annotated_by_genes_",
																											format(Sys.time(), "%y%m%d%H%M"),
																											".tsv")))
			
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
			
			req(annot_output_file()) 
			req(input$ped_tsv)
			
			pedigree_file <- isolate(input$ped_tsv$datapath)
			required_overlap <- isolate(input$th_cnv)
			inherit_output_file(file.path(results_dir, paste0("cnvs_inheritance_",
																												format(Sys.time(), "%y%m%d%H%M"),
																												".tsv")))
			
			ret <- MCNV2::compute_inheritance(cnvs_file = annot_output_file(), 
																				pedigree_file = pedigree_file, 
																				output_file = inherit_output_file(), 
																				overlap = required_overlap)
			
			incProgress(0.5, detail = "Inheritance calculation...")
			
			if(ret == 0){
				inherit_check(TRUE)
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
	
	observeEvent(inherit_output_file(), {
		if (file.exists(inherit_output_file())) {
			updateActionButton(session, "submit_mpviz", disabled = FALSE)
		} else {
			updateActionButton(session, "submit_mpviz", disabled = TRUE)
		}
	})
	
	output$conditional_input <- renderUI({
		if (inherit_check()) {
			tagList(
				tags$p(paste("✅ Loaded inheritance file:\n",
										 basename(inherit_output_file())))
			)
		} else {
			tagList(
				fileInput("preproc_file", label = "Preprocessed input (mandatory)"),
				verbatimTextOutput("preproc_tsv_status")
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
			
			internal_columns <- c("Chr","Start","End", "Exon_Overlap",
														"t_Start","t_End","STOP", "LOEUF",
														"segmentaldup_Overlap", "Size",
														"cnv_problematic_region_overlap", 
														"transcript_length")
			
			spec <- readr::spec_tsv(inherit_output_file())
			
			num_cols <- names(spec$cols)[
				vapply(spec$cols, function(x) {
					inherits(x, "collector_double") || inherits(x, "collector_number")
				}, logical(1))
			]
			num_cols <- num_cols[!num_cols %in% internal_columns]
			if(length(num_cols) > 0){
				quality_metrics(num_cols)
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
			tagList(fluidRow(
				column(6,
							 numericInput("quality_metric_min", "Min:", value = rng_infos$min, 
							 						 min = rng_infos$min, max = rng_infos$max, 
							 						 width = "100%")
				),
				column(6,
							 numericInput("quality_metric_max", "Max:", value = rng_infos$max, 
							 						 min = rng_infos$min, max = rng_infos$max, 
							 						 width = "100%")
				)
			))
			
		} else {
			tags$p("❌ The chosen quality metric is not numeric.")
		}
		
	})
	
	
	# When min changes, ensure max is >= min
	observeEvent(input$min_val, {
		if (input$max_val < input$min_val) {
			updateNumericInput(session, "max_val", value = input$min_val)
		}
		updateNumericInput(session, "max_val", min = input$min_val)
	})
	
	# When max changes, ensure min is <= max
	observeEvent(input$max_val, {
		if (input$min_val > input$max_val) {
			updateNumericInput(session, "min_val", value = input$max_val)
		}
		updateNumericInput(session, "min_val", max = input$max_val)
	})
	
	output$plots_overview <- renderUI({
		req(inherit_output_file())
		req(input$quality_metric)
		
		outputs <- list(br(),
										plotlyOutput(outputId = "plot_overview_del", 
																 height = "500px", width = "100%"), hr(), 
										plotlyOutput(outputId = "plot_overview_dup", 
																 height = "500px", width = "100%"), hr())
		
		return(outputs)
	})
	
	observeEvent(input$submit_mpviz, {
		withProgress(message = "Running analysis...", value = 0, {
			req(inherit_output_file())
			req(input$quality_metric)
			
			cnv_range <- isolate(input$cnv_range)
			# convert CNV size range
			min_cnv_size <- parse_cnv_size_value(cnv_range[1])
			max_cnv_size <- parse_cnv_size_value(cnv_range[2])
			exclusion_list <- isolate(exclusion_list())
			min_loeuf <- isolate(input$min_loeuf)
			min_transcript_ov <- isolate(input$min_transcript_ov)
			quality_metric <- isolate(input$quality_metric)
			max_prb_region_ov <- isolate(input$max_prb_region_ov)
			# get quality metric thresholds
			min_qty_metric <- isolate(input$quality_metric_min)
			max_qty_metric <- isolate(input$quality_metric_max)
			plot_type <- isolate(input$plot_type)
			
			dataset <- arrow::open_tsv_dataset(inherit_output_file())
			
			dataset <- dataset %>% mutate(Size_Range = case_when(
				Size < 30000 ~ "1-30kb",
				Size < 50000 ~ "30-50kb",
				Size < 100000 ~ "50-100kb",
				Size < 200000 ~ "100-200kb",
				Size < 500000 ~ "200-500kb",
				Size < 1000000 ~ "500-1M",
				Size >= 1000000 ~ ">1M",
				TRUE ~ NA_character_
			))
			
			# calcul de la MP globale
			transmission_col <- ifelse(test = isolate(input$transmission) == 1, 
																 yes = "transmitted_cnv", 
																 no = "transmitted_gene")
			
			if(transmission_col == "transmitted_gene"){
				# suppress intergenic CNV
				dataset <- dataset %>% 
					dplyr::filter(transmitted_gene != "intergenic")
			}
			
			incProgress(0.5, detail = "Create subset...")
			
			filtered_dataset <- dataset %>% 
				dplyr::filter(
					Size >= min_cnv_size,
					Size <= max_cnv_size
				)
			
			if(min_transcript_ov > 0){
				filtered_dataset <- filtered_dataset %>% 
					dplyr::filter(
						bp_overlap >= min_transcript_ov/100
					)
			}
			
			if(max_prb_region_ov < 100){
				filtered_dataset <- filtered_dataset %>% 
					dplyr::filter(
						cnv_problematic_region_overlap <= max_prb_region_ov/100
					)
			}
			
			# apply exclusion list
			if(!is.null(exclusion_list)){
				# get cnv_id that contain an excluded genes
				cnv_ids_2_remove <- unique(filtered_dataset %>% select(cnv_id, GeneID) %>% 
																	 	filter(GeneID %in% exclusion_list) %>% pull(cnv_id))
				
				if(length(cnv_ids_2_remove) > 0){
					filtered_dataset <- filtered_dataset %>% dplyr::filter(!cnv_id %in% cnv_ids_2_remove)
				}
			}
			
			# apply LOEUF filter
			if(min_loeuf > 0){
				# get cnv_id that contain an excluded genes (LOEUF < min_loeuf)
				cnv_ids_2_remove <- unique(filtered_dataset %>% select(cnv_id, LOEUF) %>% 
																	 	filter(LOEUF < min_loeuf) %>% pull(cnv_id))
				
				if(length(cnv_ids_2_remove) > 0){
					filtered_dataset <- filtered_dataset %>% dplyr::filter(!cnv_id %in% cnv_ids_2_remove)
				}
			}
			
			# update reactive filtered dataset
			filtered_ds(filtered_dataset)
			
			# TODO: control this calculation
			mp <- filtered_ds() %>% 
				dplyr::select(cnv_id, Type, all_of(transmission_col)) %>%
				mutate(transmission = as.logical(.data[[transmission_col]])) %>%
				distinct() %>%
				group_by(Type, transmission) %>% 
				summarise(n = n(), .groups = "drop") %>%
				collect() %>% 
				tidyr::complete(Type, transmission = c(FALSE, TRUE), fill = list(n = 0)) %>%
				group_by(Type) %>% 
				summarise(n_de_novo = sum(n[!transmission]), 
									n_inherited = sum(n[transmission]),
									prop_inherited = n_inherited / (n_de_novo + n_inherited))
			
			n_cnvs_all <- dataset %>% dplyr::select(cnv_id) %>% distinct() %>% summarize(n = n()) %>% collect()
			n_cnvs_filtered <- filtered_ds() %>% dplyr::select(cnv_id) %>% distinct() %>% summarize(n = n()) %>% collect()
			
			output$n_filtered_cnvs <- renderValueBox({
				valueBox(
					format(n_cnvs_filtered$n, big.mark = ","), 
					paste0("Filtered CNVs over ", format(n_cnvs_all$n, big.mark = ",")), 
					icon = icon("filter"), color = "purple"
				)
			})
			
			output$mp_del <- renderValueBox({
				valueBox(
					ifelse(test = !is.null(mp[mp$Type == "DEL",]$prop_inherited),
								 yes = pct(round(mp[mp$Type == "DEL",]$prop_inherited, digits = 3)),
								 no = "NA"), "Global MP (DEL)", icon = icon("square-minus"), 
					color = "red"
				)
			})
			
			output$mp_dup <- renderValueBox({
				valueBox(
					ifelse(test = !is.null(mp[mp$Type == "DUP",]$prop_inherited),
								 yes = pct(round(mp[mp$Type == "DUP",]$prop_inherited, digits = 3)),
								 no = "NA"), "Global MP (DUP)", icon = icon("square-plus"), 
					color = "teal"
				)
			})
			
			incProgress(0.7, detail = "Plotting deletions in progress...")
			output$plot_overview_del <- renderPlotly({
				
				if(plot_type == "mp_quality"){
					dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(Type == "DEL"),
																						quality_metric = quality_metric,
																						thresholds = c(min_qty_metric, max_qty_metric),
																						transmission_col = transmission_col),
											 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(Type == "DEL"),
											 						 quality_metric = quality_metric,
											 						 thresholds = c(min_qty_metric, max_qty_metric),
											 						 transmission_col = transmission_col))
				} else {
					dat <- mp_by_size(ds = filtered_ds() %>% dplyr::filter(Type == "DEL"),
														transmission_col = transmission_col)
				}
				
				
				if(nrow(dat) > 0){
					
					if(plot_type == "mp_quality"){
						p <- plot_mp_vs_metric(dt = dat,
																	 title = "Mendelian Precision - DEL",
																	 subtitle = "All filtered CNVs",
																	 y_lab = "Mendelian Precision",
																	 x_lab = paste0(quality_metric, " threshold (≥)"))
					} else {
						p <- plot_mp_vs_size(dt = dat,
																 title = "Mendelian Precision - DEL",
																 subtitle = "All filtered CNVs",
																 y_lab = "Mendelian Precision",
																 x_lab = "CNV Size (binned)")
					}
					
					baseline_DEL_plot(p)
					ggplotly(p)
				} else {
					NULL
				}
				
			})
			
			incProgress(0.9, detail = "Plotting duplications in progress...")
			output$plot_overview_dup <- renderPlotly({
				
				if(plot_type == "mp_quality"){
					dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(Type == "DUP"),
																						quality_metric = quality_metric,
																						thresholds = c(min_qty_metric, max_qty_metric),
																						transmission_col = transmission_col),
											 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(Type == "DUP"),
											 						 quality_metric = quality_metric,
											 						 thresholds = c(min_qty_metric, max_qty_metric),
											 						 transmission_col = transmission_col))
				} else {
					dat <- mp_by_size(ds = filtered_ds() %>% dplyr::filter(Type == "DUP"),
														transmission_col = transmission_col)
				}
				
				
				if(nrow(dat) > 0){
					
					if(plot_type == "mp_quality"){
						p <- plot_mp_vs_metric(dt = dat,
																	 title = "Mendelian Precision - DUP",
																	 subtitle = "All filtered CNVs",
																	 y_lab = "Mendelian Precision",
																	 x_lab = paste0(quality_metric, " threshold (≥)"))
					} else {
						p <- plot_mp_vs_size(dt = dat,
																 title = "Mendelian Precision - DUP",
																 subtitle = "All filtered CNVs",
																 y_lab = "Mendelian Precision",
																 x_lab = "CNV Size (binned)")
					}
					
					baseline_DUP_plot(p)
					ggplotly(p)
				} else {
					NULL
				}
				
			})
			
			#TODO: modify so the completion natification appears AFTER plots
			showNotification("Plotting complete!", type = "message")
			
		})
	})
	
	output$filtered_tbl <- renderDataTable({
		req(filtered_ds())
		dat <- filtered_ds() %>% head(1000) %>% collect()
		datatable(dat, options = list(scrollX = TRUE), rownames = FALSE)
	})
	
	output$ddl_filtered_tbl <- downloadHandler(
		filename = function() {
			paste0("MCNV2_filtered_data_", Sys.Date(), ".csv")
		},
		content = function(file) {
			withProgress(
				message = 'Preparing large dataset...',
				detail = 'Please wait while we generate your TSV.',
				value = 0, {
					
					incProgress(0.2, detail = "Collecting data...")
					dat <- filtered_ds() %>% collect()
					
					incProgress(0.6, detail = "Writing TSV file...")
					write.csv(dat, file, row.names = FALSE)
					
					incProgress(1, detail = "Done!")
				}
			)
		}
	)
	
	observeEvent(filtered_ds(), {
		if (!is.null(filtered_ds())) {
			updateActionButton(session, "goto_finetuning", disabled = FALSE)
		} else {
			updateActionButton(session, "goto_finetuning", disabled = TRUE)
		}
	})
	
	observeEvent(input$goto_finetuning, {
		updateTabItems(session, "tabs", "fine_tuning")  # move to the "mp_exploration" tab
	})
	
	output$finetune_ui <- renderUI({
		req(filtered_ds())
		
		outputs <- lapply(quality_metrics(), function(metric) {
			
			col_vals <- filtered_ds() %>% pull(all_of(metric))
			rng <- range(col_vals, na.rm = TRUE)
			
			ns_metric <- paste0("filter_", metric)
			
			fluidRow(
				class = "align-items-end",
				column(4,
							 pickerInput(
							 	inputId = paste0(ns_metric, "_op"),
							 	label = metric,
							 	choices = c("-" = "NA", "≥" = ">=", "≤" = "<="),
							 	selected = "NA",
							 	options = list(
							 		style = "btn-light",   # Bootstrap color theme (btn-default, btn-light, btn-info, etc.)
							 		size = "sm",             # sm, lg
							 		`dropup-auto` = FALSE    # optional: keep dropdown direction stable
							 	)
							 )
				),
				column(8,
							 numericInput(
							 	inputId = paste0(ns_metric, "_val"),
							 	label = "value",
							 	value = rng[1],
							 	min = rng[1],
							 	max = rng[2]
							 )
				)
			)
		})
		
		tagList(
			h4("Filters"),
			radioButtons("cnv_type", "CNV type", 
									 choices = c("DEL","DUP"), selected = "DEL"),
			outputs,
			hr(),
			actionButton("submit_ft_mpviz", label = "Apply filters",
									 icon = icon("gear"), disabled = FALSE)
		)
	})
	
	observe({
		lapply(quality_metrics(), function(metric) {
			ns_metric <- paste0("filter_", metric)
			op_id <- paste0(ns_metric, "_op")
			val_id <- paste0(ns_metric, "_val")
			
			observeEvent(input[[op_id]], {
				if (is.null(input[[op_id]])) return(NULL)
				if (input[[op_id]] == "NA") {
					shinyjs::disable(val_id)
				} else {
					shinyjs::enable(val_id)
				}
			}, ignoreInit = TRUE)
		})
	})
	
	# --- activer bouton seulement si au moins un filtre est actif
	observe({
		# On construit un vecteur des statuts de tous les filtres
		filter_states <- sapply(quality_metrics(), function(metric) {
			op <- input[[paste0("filter_", metric, "_op")]]
			val <- input[[paste0("filter_", metric, "_val")]]
			if (is.null(op) || op == "NA") return(FALSE)
			TRUE
		})
		
		# Si au moins un TRUE → activer le bouton
		if (any(filter_states)) {
			shinyjs::enable("submit_ft_mpviz")
		} else {
			shinyjs::disable("submit_ft_mpviz")
		}
	})

	
	observeEvent(input$submit_ft_mpviz, {
		withProgress(message = "Running analysis...", value = 0, {
			req(inherit_output_file())
			req(filtered_ds())
			cnv_type <- isolate(input$cnv_type)
			plot_type <- isolate(input$plot_type)
			quality_metric <- isolate(input$quality_metric)
			# get quality metric thresholds
			min_qty_metric <- isolate(input$quality_metric_min)
			max_qty_metric <- isolate(input$quality_metric_max)
			# transmission type
			transmission_col <- ifelse(test = isolate(input$transmission) == 1, 
																 yes = "transmitted_cnv", 
																 no = "transmitted_gene")
			
			### P1: before_add_filters ###
			output$before_add_filters <- renderPlotly({
				if(cnv_type == "DEL"){
					p <- baseline_DEL_plot()
				} else {
					p <- baseline_DUP_plot()
				}
				clean_plot_for_plotly(p)
			})
			
			### P2 - 4: after_add_filters ###
			# get filters
			df <- filtered_ds()
			
			for (m in quality_metrics()) {
				op <- isolate(input[[paste0("filter_", m, "_op")]])
				val <- isolate(input[[paste0("filter_", m, "_val")]])
				
				if (!is.null(op) && op != "NA" && !is.na(val)) {
					
					if (op == ">=") {
						df <- df %>% filter(.data[[m]] >= val)
					} else if (op == "<=") {
						df <- df %>% filter(.data[[m]] <= val)
					}
				}
			}
			
			
			output$after_add_filters <- renderPlotly({
				p <- NULL
				
				# Aucun chamgement
				if((filtered_ds() %>% count %>% collect) == (df %>% count %>% collect)){
					print("aucun changement")
					if(cnv_type == "DEL"){
						p <- baseline_DEL_plot()
					} else {
						p <- baseline_DUP_plot()
					}
					return(clean_plot_for_plotly(p))
				}
				
				# if any changes
				if(plot_type == "mp_quality"){
					dat <- rbind(mp_vs_metric_by_size(ds =  df %>% dplyr::filter(Type == cnv_type),
																						quality_metric = quality_metric,
																						thresholds = c(min_qty_metric, max_qty_metric),
																						transmission_col = transmission_col),
											 mp_vs_metric(ds = df %>% dplyr::filter(Type == cnv_type),
											 						 quality_metric = quality_metric,
											 						 thresholds = c(min_qty_metric, max_qty_metric),
											 						 transmission_col = transmission_col))
				} else {
					dat <- mp_by_size(ds = df %>% dplyr::filter(Type == cnv_type),
														transmission_col = transmission_col)
				}
				
				if(nrow(dat) > 0){
					
					if(plot_type == "mp_quality"){
						p <- plot_mp_vs_metric(dt = dat,
																	 title = paste0("Mendelian Precision - ", cnv_type),
																	 subtitle = "All filtered CNVs",
																	 y_lab = "Mendelian Precision",
																	 x_lab = paste0(quality_metric, " threshold (≥)"))
					} else {
						p <- plot_mp_vs_size(dt = dat,
																 title = paste0("Mendelian Precision - ", cnv_type),
																 subtitle = "All filtered CNVs",
																 y_lab = "Mendelian Precision",
																 x_lab = "CNV Size (binned)")
					}
					
				}
				
				clean_plot_for_plotly(p)
			})
			
		})
	})
	
	
}
