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
library(fresh)

my_theme <- create_theme(
	theme = "default", # Can also use a Bootswatch template like "cerulean"
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
		menuItem("MP Exploration", icon = icon("chart-line"), tabName = "mp_exploration"),
		menuItem("De novo Analysis", icon = icon("magnifying-glass-chart"), tabName = "denovo_analysis")
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
								h4("Parameters"),
								selectInput("build", label = "Genome build",
														choices = list("GRCh38/hg38" = 38, "GRCh37/hg19" = 37),
														selected = 38),
								#radioButtons("build", "Genome build", c("38","37"), selected = "38", inline = TRUE),
								numericInput("th_prob", "Problematic regions threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
								numericInput("th_cnv", "Inheritance threshold (child CNV proportion)", 0.50, min = 0, max = 1, step = 0.05),
								hr(),
								div(class = "run-status", uiOutput("run_state")),
								div(class = "outdir-box", verbatimTextOutput("outdir_display")),
								actionButton("submit_preprocess", label = "Submit",
														 icon = icon("gear"), disabled = TRUE)
								
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
																 ))               
							)
						)
		),
		
		tabItem(tabName = "mp_exploration",
						h3("MP Exploration"),
						helpText(
							"This step requires a preprocessed file created at the ",
							tags$a("Preprocessing step", onclick = "openTab('prepocessing')", href = "#")
						),
						br(),
						sidebarLayout(
							sidebarPanel(
								h4("Input"),
								uiOutput("conditional_input"),
								hr(),
								h4("Parameters"),
								radioButtons("transmission", label = "Transmission type", 
														 choices = list("CNV level" = 1, "Gene level" = 2), 
														 selected = 2, inline = TRUE),
								checkboxGroupInput("cnv_type", label = "CNV type", 
														choices = c("DEL", "DUP"), 
														selected = c("DEL","DUP"), inline = TRUE),
								selectizeInput("quality_metric", label = "Quality metric",
															 options = NULL, choices = NULL),
								hr(),
								h4("Inclusion criteria"),
								uiOutput("score_ui"),
								uiOutput("size_ui"),
								sliderInput("min_exon_ov", "% Exon overlap minimum",
														min = 0, max = 100, value = 0, step = 10),
								shinyWidgets::sliderTextInput(inputId = "cnv_size", 
																							label = "CNV size range", 
																							choices = c(1,5,10,15,20,25,30)),
								hr(),
								h4("Exlusion criteria"),
								uiOutput("score_ui"),
								hr(),
								actionButton("submit_display", label = "Apply filters",
														 icon = icon("gear"), disabled = TRUE)
								
							),
							
							# Show a plot of the generated distribution
							mainPanel(
								tableOutput("summary_table"),
								tableOutput("header"),
							)
						)
		),
		
		tabItem(tabName = "denovo_analysis",
						h3("De Novo Analysis"),
						helpText("HELP"),
						br()
		)
	)
)

dashboardPage(skin = "black",
							dashboardHeader(title = "MCNV2 - Mendelian CNV Validation"),
							sidebar,
							body)
