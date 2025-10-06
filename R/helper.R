#' @export
check_input_file <- function(filepath, required_cols = NULL, sep = "\t") {
	
	# Vérifie que le fichier existe
	if (!file.exists(filepath)) {
		return(list(status = FALSE, msg = paste0("❌ File not found: ", filepath)))
	}
	
	# Essaie de lire juste l'en-tête (pour aller vite)
	header <- names(read.table(filepath, header = TRUE, sep = sep, nrows = 1, check.names = FALSE))
	
	# Si des colonnes sont requises, on les teste
	if (!is.null(required_cols)) {
		missing_cols <- setdiff(required_cols, header)
		if (length(missing_cols) > 0) {
			return(list(status = FALSE, msg = paste0("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))))
		}
	}
	
	return(list(status = TRUE, msg = "✅ File check passed"))
	
}

