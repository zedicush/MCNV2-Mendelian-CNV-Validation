#' @export
setup_python_env <- function(envname = "r-MCNV2") {
	library(reticulate)
	
	if (!virtualenv_exists(envname)) virtualenv_create(envname)
	
	req_file <- system.file("python/requirements.txt", package = "MCNV2")
	cmd <- paste(shQuote(path.expand(file.path("~", ".virtualenvs", envname, "bin", "python"))),
							 "-m pip install --upgrade -r", shQuote(req_file))
	system(cmd, intern = TRUE)
	
	# Tell reticulate to use this env
	use_virtualenv(envname, required = TRUE)
	invisible(TRUE)
}
