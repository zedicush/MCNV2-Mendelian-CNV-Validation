# R/run_app.R
#' Run the Shiny Application
#'
#' @export
launch <- function() {
	appDir <- system.file("shiny-app", package = "MCNV2")
	if (appDir == "") {
		stop("Could not find shiny-app directory. Try re-installing `MCNV2`.", call. = FALSE)
	}
	shiny::runApp(appDir, display.mode = "normal",
								launch.browser = TRUE)
}
