#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

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
	
	prob_regions_file <- reactiveVal(NULL)
	annot_output_file <- reactiveVal(NULL)
	inherit_output_file <- reactiveVal(NULL)
	
	filtered_ds <- reactiveVal(NULL) #arrow dataset to avoid loading file in memory
	
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
			
			pedigree_file <- isolate(input$ped_tsv$datapath)
			required_overlap <- isolate(input$th_cnv)
			inherit_output_file(file.path(results_dir, "cnvs_inheritance.tsv"))
			
			ret <- MCNV2::compute_inheritance(cnvs_file = annot_output_file(), 
																				pedigree_file = pedigree_file, 
																				output_file = inherit_output_file(), 
																				overlap = required_overlap)
			
			incProgress(0.5, detail = "Inheritance calculation...")
			
			if(ret == 0){
				inherit_check(TRUE)
				#inherit_output_file("~/workspace/projects/Support/BIOINF-219/data/DEL_DUP_inheritance_readyForShiny_annot.tsv")
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
														"t_Start","t_End","STOP", 
														"LOEUF", "segmentaldup_Overlap")
			
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
																 height = "500px"), hr(), 
										plotlyOutput(outputId = "plot_overview_dup", 
																 height = "500px"), hr())
		
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
			min_exon_ov <- isolate(input$min_exon_ov)
			min_transcript_ov <- isolate(input$min_transcript_ov)
			exclus_genes <- isolate(input$exclus_genes)
			quality_metric <- isolate(input$quality_metric)
			# get quality metric thresholds
			min_qty_metric <- isolate(input$quality_metric_min)
			max_qty_metric <- isolate(input$quality_metric_max)
			
			dataset <- arrow::open_tsv_dataset(inherit_output_file())
			
			dataset <- dataset %>% mutate(Size_Range = case_when(
				Size < 30000 ~ "1-30kb",
				Size < 40000 ~ "30-40kb",
				Size < 50000 ~ "40-50kb",
				Size < 100000 ~ "50-100kb",
				Size < 200000 ~ "100-200kb",
				Size < 500000 ~ "200-500kb",
				Size < 1000000 ~ "500-1M",
				Size >= 1000000 ~ ">1M",
				TRUE ~ NA_character_
			))
			
			# calcul de la MP globale
			transmission_col <- ifelse(test = input$transmission == 1, 
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
					#Exon_Overlap >= min_exon_ov/100,
					bp_overlap >= min_transcript_ov/100,
					Size >= min_cnv_size,
					Size <= max_cnv_size,
					.data[[input$quality_metric]] >= min_qty_metric,
					.data[[input$quality_metric]] <= max_qty_metric
				)
			
			# apply exclusion list
			# TODO: est-ce qu'on retire un CNV dès qu'il impact un "green" genes?
			if(!is.null(exclus_genes)){
				lines <- tryCatch(readLines(exclus_genes$datapath, warn = FALSE), 
													error = function(e) character())
				exclusion_list <- unique(trimws(lines[nzchar(trimws(lines))]))
				filtered_dataset <- filtered_dataset %>% dplyr::filter(!GeneID %in% exclusion_list)
			}		
			
			# update reactive filtered dataset
			filtered_ds(filtered_dataset)

			
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
				dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(Type == "DEL"), 
																					quality_metric = quality_metric, 
																					transmission_col = transmission_col), 
										 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(Type == "DEL"), 
										 						 quality_metric = quality_metric, 
										 						 transmission_col = transmission_col))
				if(nrow(dat) > 0){
					p <- plot_mp_vs_metric(dt = dat, 
																 title = "Mendelian Precision - DEL", 
																 subtitle = "All filtered CNVs", 
																 y_lab = "Mendelian Precision", 
																 x_lab = paste0(quality_metric, " threshold (≥)"))
					
					ggplotly(p) 
				} else {
					NULL
				}
				
			}) 
			
			incProgress(0.9, detail = "Plotting duplications in progress...")
			output$plot_overview_dup <- renderPlotly({ 
				dat <- rbind(mp_vs_metric_by_size(ds =  filtered_ds() %>% dplyr::filter(Type == "DUP"), 
																					quality_metric = quality_metric, 
																					transmission_col = transmission_col), 
										 mp_vs_metric(ds = filtered_ds() %>% dplyr::filter(Type == "DUP"), 
										 						 quality_metric = quality_metric, 
										 						 transmission_col = transmission_col))
				
				if(nrow(dat) > 0){
					p <- plot_mp_vs_metric(dt = dat, 
																 title = "Mendelian Precision - DUP", 
																 subtitle = "All filtered CNVs", 
																 y_lab = "Mendelian Precision", 
																 x_lab = paste0(quality_metric, " threshold (≥)"))
					
					ggplotly(p) 
				} else {
					NULL
				}
			}) 
			showNotification("Plotting complete!", type = "message")
			
		})
	})
	
	output$investigate_del <- renderUI({
		req(inherit_output_file())
		req(input$quality_metric)
		
		plot_types <- c(
			
		)
		outputs <- list()
		
		
		
		return(outputs)
	})
	
}
