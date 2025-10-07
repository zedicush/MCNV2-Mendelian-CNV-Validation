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
	inheritance_output_file <- reactiveVal(NULL)
	
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
				inheritance_output_file("~/workspace/projects/Support/BIOINF-219/data/DEL_DUP_inheritance_readyForShiny.tsv")
				output$inherit_tsv_status <- renderText(paste("✅ Path to inheritance file:\n",
																											inheritance_output_file()))
				dat <- readr::read_tsv(inheritance_output_file(), n_max = 50, show_col_types = FALSE)
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
				tags$text(paste("✅ Path to inheritance file:\n",
												inheritance_output_file()))
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
			inheritance_output_file(filepath)
			output$preproc_tsv_status <- renderText(ret$msg)
		} else {
			output$preproc_tsv_status <- renderText(ret$msg)
		}
		
	})
	
	observeEvent(inheritance_output_file(), {
		if (file.exists(inheritance_output_file())) {
			dat <- readr::read_tsv(inheritance_output_file(), n_max = 1, 
														 show_col_types = FALSE)
			updateSelectizeInput(inputId = "quality_metric", 
													 choices = colnames(dat))
		} else {
			updateSelectizeInput(inputId = "quality_metric", choices = NULL)
		}
	})
	
	
	
}
