cnv_size_colors <- c(
	"1-30kb" = "#A6CEE3",
	"30-50kb" = "#1F78B4",
	"50-100kb" = "#33A02C",
	"100-200kb" = "#6A3D9A",
	"200-500kb" = "#CAB2D6",
	"500-1M" = "#FDBF6F",
	">1M" = "#FF7F00",
	"All" = "#E31A1C"
)

cnv_size_shapes <- c(
	"1-30kb" = "1",
	"30-50kb" = "2",
	"50-100kb" = "3",
	"100-200kb" = "4",
	"200-500kb" = "5",
	"500-1M" = "6",
	">1M" = "7",
	"All" = "8"
)

#' @export
check_input_file <- function(filepath, 
														 file_type = c("cnv", "ped","prb", "preproc"), 
														 sep = "\t") {
	file_type <- match.arg(file_type)
	
	if (!file.exists(filepath)) {
		return(list(status = FALSE, msg = paste0("❌ File not found: ", filepath)))
	}
	
	# Read header
	header <- names(read.table(filepath, header = TRUE, sep = sep, 
														 nrows = 1, check.names = FALSE))
	header_lower <- tolower(header)
	
	# Define accepted column variants per file type
	schema_list <- list(
		cnv = list(
			ordered = list(
				chr   = c("chr", "chrom", "seqnames"),
				start = c("start", "begin", "pos_start"),
				end   = c("end", "stop", "pos_end")
			),
			unordered = list(
				sample = c("sampleid", "sample_id", "id_sample"),
				type   = c("type", "svtype", "variant_type")
			)
		),
		ped = list(
			ordered = list(
				iid = c("iid", "individual", "sampleid", "childid", "id")
			),
			unordered = list(
				father = c("pid", "father", "dad", "fatherid"),
				mother = c("mid", "mother", "mom", "motherid")
			)
		),
		prb = list(
			ordered = list(
				chr   = c("chr", "chrom", "seqnames"),
				start = c("start", "begin", "pos_start"),
				end   = c("end", "stop", "pos_end")
			),
			unordered = list()
		),
		preproc = list(
			ordered = list(),
			unordered = list(
				cnvid   = c("cnv_id"),
				chr   = c("chr", "chrom", "seqnames","chr_child", "chrom_child", "seqnames_child"),
				start = c("start", "begin", "pos_start","start_child", "begin_child", "pos_start_child"),
				end   = c("end", "stop", "pos_end","end_child", "stop_child", "pos_end_child"),
				sample = c("sampleid", "sample_id", "id_sample","sampleid_child"),
				type   = c("type", "svtype", "variant_type","type_child"),
				length = c("size","length","size_child"),
				transcript_overlap = c("transcript_overlap","transcript_bp_overlap","bp_overlap"),
				geneid = c("gene_id","geneid"),
				biotype = c("biotype","gene_type"),
				cnv_inherit = c("inheritance_by_cnv","transmitted_cnv"),
				gene_inherit = c("inheritance_by_gene","transmitted_gene"),
				loeuf = c("loeuf")
			)
		),
		other = list(
			ordered = list(),
			unordered = list()
		)
	)
	
	schema <- schema_list[[file_type]]
	
	# --- Check ordered columns ---
	if (length(schema$ordered) > 0) {
		required_n <- length(schema$ordered)
		if (length(header_lower) < required_n) {
			return(list(status = FALSE, msg = paste0("❌ File has fewer than ", required_n, " columns.")))
		}
		
		ordered_names <- names(schema$ordered)
		for (i in seq_along(ordered_names)) {
			colname <- ordered_names[i]
			variants <- schema$ordered[[colname]]
			if (!header_lower[i] %in% variants) {
				return(list(
					status = FALSE,
					msg = paste0("❌ Column ", i, " must be one of: ", paste(variants, collapse = ", "))
				))
			}
		}
	}
	
	# --- Check unordered columns ---
	if (length(schema$unordered) > 0) {
		missing_cols <- c()
		for (colname in names(schema$unordered)) {
			variants <- schema$unordered[[colname]]
			if (!any(header_lower %in% variants)) {
				missing_cols <- c(missing_cols, colname)
			}
		}
		if (length(missing_cols) > 0) {
			return(list(status = FALSE, msg = paste0("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))))
		}
	}
	
	return(list(status = TRUE, msg = "✅ File check passed!"))
}

#' @export
load_gene_ids <- function(filepath) {
	
	# --- Step 1: basic existence check
	if (!file.exists(filepath)) {
		return(list(status = FALSE, 
								msg = paste0("❌ File not found: ", filepath),
								genes = NULL))
	}
	
	# --- Step 2: read safely, only first column, ignore empty lines or comments
	df <- NULL
	df <- tryCatch(
		readr::read_tsv(
			filepath,
			col_names = FALSE,
			col_types = readr::cols(.default = readr::col_character()),
			comment = "#",
			progress = FALSE,
			show_col_types = FALSE
		),
		error = function(e) print(e))
	
	if(is.null(df)) return(list(status = FALSE, 
															msg = paste0("❌ Unable to read: ", filepath),
															genes = NULL))
															
	# --- Step 3: structure checks
	if (ncol(df) < 1) return(list(status = FALSE, 
																msg = "❌ File must contain at least one column.",
																genes = NULL))
	if (nrow(df) == 0) return(list(status = FALSE, 
																 msg = "❌ File is empty.",
																 genes = NULL))

	# --- Step 4: keep only first column
	genes <- df[[1]]
	
	# --- Step 5: validate format (Ensembl Gene IDs)
	valid_pattern <- "^ENSG[0-9]+$"
	invalid <- genes[!grepl(valid_pattern, genes)]
	if (length(invalid) > 0) {
		return(list(status = FALSE, 
								msg = "❌ Invalid Ensembl Gene IDs detected.",
								genes = NULL))
	}
	
	# --- Step 6: return clean vector (no duplicates, no NAs)
	genes <- unique(na.omit(trimws(genes)))

	return(list(status = TRUE, 
							msg = paste0("✅ Loaded ", length(genes), " IDs."),
							genes = genes))
}


#' @export
annotate <- function(cnvs_file, prob_regions_file, output_file, 
										 genome_version = 38, bedtools_path){
	
	gene_annotation_script <- system.file("python",
																				"gene_annotation.py", 
																				package = "MCNV2")
	gene_resource_file <- system.file("resources",
																		"gene_resources.tsv", 
																		package = "MCNV2")
	
	args <- c(
		gene_annotation_script,
		"--cnv", cnvs_file,
		"--gene_resource", gene_resource_file,
		"--prob_regions", prob_regions_file,
		"--out", output_file,
		"--genome_version", genome_version,
		"--bedtools_path", bedtools_path
	)
	
	ret <- system2(
		command = reticulate::py_exe(),  # ex: ~/.virtualenvs/r-MCNV2/bin/python
		args = args, 
		stdout = TRUE,  # capture stdout si tu veux l'afficher
		stderr = TRUE
	)
	
	if(!file.exists(output_file)){
		return(1)
	} else {
		return(0)
	}
}

#' @export
compute_inheritance <- function(cnvs_file, pedigree_file, output_file, 
																overlap = 0.5){

	compute_inheritance_script <- system.file("python", "compute_inheritance.py", 
																				package = "MCNV2")
	
	args <- c(
		compute_inheritance_script, 
		"--cnv_geneAnnot", cnvs_file,
		"--pedigree", pedigree_file, 
		"--output", output_file, 
		"--overlap", overlap
	)
	
	ret <- system2(
		command = reticulate::py_exe(),  # ex: ~/.virtualenvs/r-MCNV2/bin/python
		args = args, 
		stdout = TRUE, 
		stderr = TRUE
	)
	
	if(!file.exists(output_file)){
		return(1)
	} else {
		return(0)
	}
}

parse_cnv_size_value <- function(x) {
	if(grepl(x = x, pattern = ">")){
		x <- 1000000000
	} else {
		x <- gsub("kb", "000", x)
		x <- gsub("Mb", "000000", x)
		x <- as.numeric(x)
	}
	
	return(x)
}

create_dynamic_ticks <- function(range_values, n_steps = 100) {
	if (!is.numeric(range_values) || length(range_values) != 2) {
		stop("range_values must be a numeric vector of length 2 (min, max).")
	}
	
	min_val <- min(range_values)
	max_val <- max(range_values)
	
	if (max_val <= 1) {
		step <- 0.1
		ticks <- seq(0, 1, by = step)
		return(list(ticks = ticks, step = step, min = min(ticks), max = max(ticks)))
	}
	
	diff_val <- max_val - min_val
	raw_step <- diff_val / n_steps
	
	# base magnitude and candidate "nice" steps (1,2,5,10 * 10^k)
	magnitude <- 10 ^ floor(log10(raw_step))
	candidates <- magnitude * c(1, 2, 5, 10)
	step <- candidates[which(candidates >= raw_step)[1]]
	if (is.na(step)) step <- tail(candidates, 1)  # fallback
	
	# force a sensible minimal step for large ranges to avoid tiny non-round steps:
	if (diff_val >= 100) {
		diff_mag <- floor(log10(diff_val))
		min_step <- 10 ^ max(0, diff_mag - 1)  # e.g. diff 192 -> min_step = 10
		if (step < min_step) step <- min_step
	}
	
	# round step to integer when >= 1
	if (step >= 1) step <- round(step)
	
	# internal round ticks (multiples of step) inside the interval
	internal_start <- ceiling(min_val / step) * step
	internal_end   <- floor(max_val / step) * step
	
	if (internal_start <= internal_end) {
		mid_ticks <- seq(internal_start, internal_end, by = step)
	} else {
		mid_ticks <- numeric(0)
	}
	
	ticks <- unique(c(min_val, mid_ticks, max_val))
	
	# tidy rounding: if step < 1 keep decimals, else integers
	digits <- if (step < 1) abs(floor(log10(step))) else 0
	ticks <- round(ticks, digits = digits)
	
	return(list(ticks = ticks, step = step, min = min(ticks), max = max(ticks)))
}

# plot_one_view <- function(df_view, type_label, mode, subtitle_tag) {
# 	if (!nrow(df_view) || !all(c("type","trans","Size_Range") %in% names(df_view))) {
# 		return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
# 	}
# 	if (identical(mode, "mp_score")) {
# 		ths <- thresholds_seq()
# 		ag <- mp_vs_score_by_size(df_view, type_label, ths)
# 		if (!nrow(ag)) return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
# 		ag$Size_Range <- factor(ag$Size_Range,
# 														levels = c("[1kb–10kb[","[10kb–30kb[","[30kb–50kb[","[50kb–100kb[",
# 																			 "[100kb–200kb[","[200kb–500kb[","[500kb–1000kb[",">=1Mb"))
# 		ggplot(ag, aes(x = threshold, y = precision, color = Size_Range)) +
# 			geom_line(linewidth = 0.9) +
# 			geom_point(size = 1.8) +
# 			geom_text(aes(label = n), vjust = -0.8, size = 2.6, alpha = 0.8, check_overlap = TRUE) +
# 			scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
# 			labs(title = paste("Mendelian Precision —", type_label),
# 					 subtitle = subtitle_tag,
# 					 x = "Score threshold (≥)", y = "Mendelian Precision") +
# 			theme_minimal(base_size = 12) +
# 			theme(plot.title = element_text(size = 13))
# 	} else {
# 		ag <- mp_by_size(df_view, type_label)
# 		if (!nrow(ag)) return(ggplot() + theme_minimal() + ggtitle(paste(type_label, "—", subtitle_tag, "(no data)")))
# 		ag$Size_Range <- factor(ag$Size_Range,
# 														levels = c("[1kb–10kb[","[10kb–30kb[","[30kb–50kb[","[50kb–100kb[",
# 																			 "[100kb–200kb[","[200kb–500kb[","[500kb–1000kb[",">=1Mb"))
# 		ggplot(ag, aes(x = Size_Range, y = mendelian_precision)) +
# 			geom_col() +
# 			geom_text(aes(label = paste0("n = ", total_CNVs))) +
# 			scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
# 			labs(title = paste("Mendelian Precision —", type_label),
# 					 subtitle = subtitle_tag,
# 					 x = "Size_Range", y = "Mendelian Precision") +
# 			theme_minimal(base_size = 8) +
# 			theme(axis.text.x = element_text(angle = 25, hjust = 1))
# 	}
# }
# 

mp_vs_metric_by_size <- function(ds, quality_metric, thresholds, 
																 transmission_col, inheritance_flag = "True") {
	ds <- isolate(ds) 
	rng_infos <- create_dynamic_ticks(thresholds, n_steps = 20)

	res <- parallel::mclapply(X = rng_infos$ticks, FUN = function(th) {
		tmp <- ds %>% dplyr::filter(.data[[quality_metric]] >= th) 
		n <- tmp %>% summarise(n()) %>% pull()
		
		if (n == 0) return(NULL)
		tmp %>%
			mutate(trans = !!sym(transmission_col) == inheritance_flag) %>%
			dplyr::group_by(Size_Range) %>%
			dplyr::summarise(n = dplyr::n(),
											 MP = round(mean(trans, na.rm = TRUE), digits = 2),
											 .groups = "drop") %>%
			dplyr::mutate(threshold = th) %>% collect()
	}, mc.cores = parallel::detectCores() - 1)
	
	dplyr::bind_rows(res)
}

mp_vs_metric <- function(ds, quality_metric, thresholds, transmission_col, 
												 inheritance_flag = "True") {
	
	ds <- isolate(ds)
	rng_infos <- create_dynamic_ticks(thresholds, n_steps = 20)
	
	res <- parallel::mclapply(X = rng_infos$ticks, FUN = function(th) {
		tmp <- ds %>% dplyr::filter(.data[[quality_metric]] >= th) 
		n <- tmp %>% summarise(n()) %>% pull()
		if (n == 0) return(NULL)
		tmp %>%
			mutate(trans = !!sym(transmission_col) == inheritance_flag) %>%
			dplyr::summarise(n = dplyr::n(),
											 MP = round(mean(trans, na.rm = TRUE), digits = 2)) %>%
			dplyr::mutate(Size_Range = "All", threshold = th) %>% collect()
	}, mc.cores = parallel::detectCores() - 1)
	
	dplyr::bind_rows(res)
}

mp_by_size <- function(ds, transmission_col, inheritance_flag = "True") {
	ds <- isolate(ds)
	
	ds %>%
		dplyr::group_by(Size_Range) %>%
		mutate(trans = !!sym(transmission_col) == inheritance_flag) %>%
		dplyr::summarise(
			n = dplyr::n(),
			MP = round(mean(trans, na.rm = TRUE), digits = 2),
			.groups = "drop"
		) %>% collect
}

plot_mp_vs_metric <- function(dt, title, subtitle, y_lab, x_lab){
	dt$Size_Range <- factor(dt$Size_Range,
													levels = c("1-30kb","30-50kb","50-100kb",
																		 "100-200kb","200-500kb","500-1M",">1M","All"))
	
	p <- ggplot(dt, aes(x = threshold, y = MP, color = Size_Range)) +
		geom_line() +
		geom_point(aes(size = n, shape = Size_Range)) +
		scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
		labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
		theme_minimal(base_size = 12) +
		scale_color_manual(values = cnv_size_colors, name = "CNV Length") +
		scale_shape_manual(values = cnv_size_shapes, name = "CNV Length")
	
	return(p)
}

plot_mp_vs_size <- function(dt, title, subtitle, y_lab, x_lab){
	dt$Size_Range <- factor(dt$Size_Range,
													levels = c("1-30kb","30-50kb","50-100kb",
																		 "100-200kb","200-500kb","500-1M",">1M","All"))
	
	p <- ggplot(dt, aes(x = Size_Range, y = MP, fill = Size_Range)) +
		geom_col() +
		labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
		geom_text(aes(label = n)) +
		theme_minimal(base_size = 12) +
		scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
		scale_fill_manual(values = cnv_size_colors, name = "CNV Length")

	return(p)
}

clean_plot_for_plotly <- function(p, 
																	x_angle = 45,      # angle pour l'axe x
																	text_size = 9,     # taille du texte des axes et geom_text
																	title_size = 13,   # taille du titre
																	legend_text_size = 9,
																	legend_title_size = 10,
																	geom_text_size = 3,
																	margin_bottom = 100) {
	
	# Redéfinit la taille par défaut des geom_text
	update_geom_defaults("text", list(size = geom_text_size))
	
	# Modifie le thème du ggplot
	p_clean <- p +
		theme_minimal(base_size = text_size) +
		theme(
			axis.text.x = element_text(
				angle = x_angle,
				hjust = ifelse(x_angle == 0, 0.5, 1),
				size = text_size
			),
			axis.text.y = element_text(size = text_size),
			plot.title = element_text(size = title_size, face = "bold"),
			legend.text = element_text(size = legend_text_size),
			legend.title = element_text(size = legend_title_size)
		)
	
	# Convertit en plotly avec marge ajustée
	ggplotly(p_clean) %>%
		layout(
			margin = list(b = margin_bottom)
		)
}


