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
library(shinyjs)
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
								fileInput("cnvs_file", label = "CNV File input (mandatory)"),
								fileInput("pedigree_file", label = "Pedigree File input (mandatory)"),
								fileInput("annot_file", label = "Annotation File input (optional)"),
								selectInput("select", label = "Genome build",
														choices = list("hg38" = 38, "hg19" = 19),
														selected = 38),
								hr(),
								actionButton("submit_preprocess", label = "Submit",
														 icon = icon("gear"), disabled = TRUE)

							),

							# Show a plot of the generated distribution
							mainPanel(
								tableOutput("summary_table"),
								tableOutput("header"),
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
								fileInput("preprocessed_file", label = "Preprocessed input (mandatory)"),
								selectInput("select", label = "Genome build",
														choices = list("hg38" = 38, "hg19" = 19),
														selected = 38),
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
