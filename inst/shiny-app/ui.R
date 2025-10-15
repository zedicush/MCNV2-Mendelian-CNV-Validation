#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(shinybusy)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)
library(fresh)
library(plotly)
library(DT)

my_theme <- create_theme(
	theme = "cosmo", # Can also use a Bootswatch template like "cerulean"
	bs_vars_button(
		default_color = "#FFF",
		default_bg = "#112446",
		default_border = "#112446",
		border_radius_base = "15px"
	),
	bs_vars_wells(
		bg = "#FFF",
		border = "#112446"
	)
)


sidebar <- dashboardSidebar(
	sidebarMenu(
		id = "tabs", 
		menuItem("Preprocessing", tabName = "prepocessing", icon = icon("users-gear")),
		menuItem("MP Exploration", tabName = "mp_exploration", icon = icon("chart-line")),
		menuItem("Fine-tuning", tabName = "fine_tuning", icon = icon("filter"))	 
	)
)

body <- dashboardBody(
	useShinyjs(), # Initialize shinyjs
	use_theme(my_theme),
	tags$script(HTML("
        var openTab = function(tabName){
          $('a', $('.sidebar')).each(function() {
            if(this.getAttribute('data-value') == tabName) {
              this.click()
            };
          });
        }
      ")),
	tags$style(type = "text/css", "
      .inline label {
        display: table-cell;
        text-align: left;
        vertical-align: middle;
      }
      .inline .form-group {
        display: table-row;
      }
    "),
	# Add a spinner in an application each time the server take more 100 milliseconds to respond.
	add_busy_spinner(spin = "fading-circle", 
									 position = "top-right", 
									 margins = c(10,10)),
	tabItems(
		tabItem(tabName = "prepocessing",
						h3("Preprocessing"),
						helpText(
							"XXXXXXXX. If you have already processed your CNV file, you can go to the ",
							tags$a("MP Exploration", onclick = "openTab('mp_exploration')", href = "#")
						),
						br(),
						sidebarLayout(
							sidebarPanel( 
								h4("Input files"),
								fileInput("cnv_tsv", label = "CNV File (mandatory, BED)"),
								verbatimTextOutput("cnv_tsv_status"),
								fileInput("ped_tsv", label = "Pedigree File (mandatory, tsv)"),
								verbatimTextOutput("ped_tsv_status"),
								fileInput("prb_tsv", label = "Problematic regions (optional, BED)"),
								verbatimTextOutput("prb_tsv_status"),
								hr(),
								h4("Parameters for inheritance calculation"),
								selectInput("build", label = "Genome build",
														choices = list("GRCh38/hg38" = 38, "GRCh37/hg19" = 37),
														selected = 38),
								numericInput("th_cnv", "Inheritance threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
								hr(),
								div(class = "run-status", uiOutput("run_state")),
								div(class = "outdir-box", verbatimTextOutput("outdir_display")),
								actionButton("submit_preprocess", label = "Submit",
														 icon = icon("gear"), disabled = TRUE),
								width = 3
							),
							mainPanel(
								conditionalPanel(condition = "input.submit_preprocess > 0",
																 bsCollapse(id = "preprocess_panel",
																 					 open = c("Annotation table (Preview)","Inheritance table (Preview)"),
																 					 multiple = TRUE,
																 					 bsCollapsePanel("Annotation table (Preview)", 
																 					 								DTOutput("preview_preproc_tbl"),
																 					 								verbatimTextOutput("annot_tsv_status"),
																 					 								hr(),
																 					 								actionButton("submit_inheritance", 
																 					 														 label = "Proceed to Inheritance calculation",
																 					 														 icon = icon("gear"), 
																 					 														 disabled = TRUE), 
																 					 								style = "info"),
																 					 bsCollapsePanel("Inheritance table (Preview)", 
																 					 								DTOutput("preview_inherit_tbl"),
																 					 								verbatimTextOutput("inherit_tsv_status"),
																 					 								hr(),
																 					 								actionButton("goto_mpexploration", 
																 					 														 label = "Go to Mendelian Precision analysis",
																 					 														 icon = icon("arrow-right"), 
																 					 														 disabled = TRUE), 
																 					 								style = "success")
																 )),
								width = 9              
							)
						)
		),
		
		tabItem(tabName = "mp_exploration",
						h3("MP Exploration"),
						helpText(
							"This step requires a preprocessed file created at the ",
							tags$a("Preprocessing step", onclick = "openTab('prepocessing')",
										 href = "#")
						),
						br(),
						sidebarLayout(
							sidebarPanel(
								h4("Input"),
								uiOutput("conditional_input"),
								hr(),
								radioButtons("transmission", label = tags$h4("Transmission type"), 
														 choices = list("CNV level" = 1, "Gene level" = 2), 
														 selected = 1, inline = TRUE),
								hr(),
								h4("CNV-level inclusion criteria"),
								helpText("Only CNVs passing filters will be used for the MP calculation."),
								sliderTextInput(
									inputId = "cnv_range",
									label = "CNV size (bp):",
									choices = c("1","30kb","40kb","50kb","100kb","200kb","500kb","1Mb",">1Mb"),
									selected = c("1", ">1Mb"),
									grid = TRUE
								),
								sliderInput("min_transcript_ov", "Min. % transcript overlap",
														min = 0, max = 100, value = 0, step = 5),
								sliderInput("max_prb_region_ov", "Max. % problematic regions overlap",
														min = 0, max = 100, value = 100, step = 5),
								hr(),
								h4("Gene-level exlusion criteria"),
								helpText("CNVs containing excluded genes (based on genes list and/or LOEUF threshold) will be excluded from the MP calculation."),
								fileInput("exclus_genes", label = "Exclusion list (Ensembl Gene IDs)"),
								verbatimTextOutput("exclus_genes_status"),
								sliderInput("min_loeuf", "Exclude genes with LOEUF <",
														min = 0, max = 1, value = 0, step = .1),
								hr(),
								h4("MP representation"),
								radioButtons("plot_type", "Plot type",
														 choices = c("MP x Quality metric" = "mp_quality",
														 						"MP x CNV size" = "mp_size"),
														 selected = "mp_quality"
								),
								conditionalPanel(condition = "input.plot_type == 'mp_quality'",
																 helpText("Choose a quantitative variable to use as a variable threshold."),
																 helpText("The MP will be calculated for CNV with a quality >= to the threshold"),
																 selectizeInput("quality_metric", label = "Quality metric",
																 							 options = NULL, choices = NULL),
																 uiOutput("qty_metric_range_ui")
								),
								actionButton("submit_mpviz", label = "Apply filters",
														 icon = icon("gear"), disabled = TRUE),
								width = 3
								
							),
							mainPanel(
								conditionalPanel(condition = "input.submit_mpviz > 0",
																 valueBoxOutput("n_filtered_cnvs"),
																 valueBoxOutput("mp_del"),
																 valueBoxOutput("mp_dup"),
																 tabsetPanel(id = "tabs",
																 						tabPanel("Overview",
																 										 uiOutput("plots_overview")
																 						),
																 						tabPanel("Filtered table",
																 										 h4("Preview"),
																 										 DTOutput("filtered_tbl"),
																 										 downloadButton("ddl_filtered_tbl", 
																 										 							 "Download CSV")
																 						)
																 ),
																 actionButton("goto_finetuning", 
																 						 label = "Go to Fine-tuning analysis",
																 						 icon = icon("arrow-right"), 
																 						 disabled = TRUE)
								),
								width = 9
							)
						)
		),
		tabItem(tabName = "fine_tuning",
						h3("Type-specific fine tuning"),
						helpText("Fine-tune your filters to maximize the MP curve"),
						br(),
						sidebarLayout(
							sidebarPanel( 
								uiOutput("finetune_ui"),
								width = 3
							),
							mainPanel(
								conditionalPanel(condition = "input.submit_ft_mpviz > 0",
																 fluidRow(
																 	column(6,
																 				 h5(tags$b("Before additional filters")),
																 				 plotlyOutput("before_add_filters")
																 	),
																 	column(6,
																 				 h5(tags$b("After additional filters")),
																 				 plotlyOutput("after_add_filters")
																 	)
																 ),
																 hr(),
																 fluidRow(
																 	column(6,
																 				 selectizeInput("p3_type", 
																 				 							 label = "Comparison 1",
																 				 							 choices = c("type1","type2","type3","type4"),
																 				 							 selected = "type1"),
																 				 plotlyOutput("p3")
																 	),
																 	column(6,
																 				 selectizeInput("p4_type", 
																 				 							 label = "Comparison 2",
																 				 							 choices = c("type1","type2","type3","type4"),
																 				 							 selected = "type2"),
																 				 plotlyOutput("p4")
																 	)
																 )
																 
								),
								width = 9              
							)
						)
		)
	)
)

dashboardPage(skin = "black",
							dashboardHeader(title = "MCNV2"),
							sidebar,
							body)
