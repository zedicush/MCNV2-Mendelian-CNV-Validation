#' @export
check_input_file <- function(filepath, 
														 file_type = c("cnv", "ped", "annot","prb", "preproc"), 
														 sep = "\t") {
	file_type <- match.arg(file_type)
	
	if (!file.exists(filepath)) {
		return(list(status = FALSE, msg = paste0("❌ File not found: ", filepath)))
	}
	
	# Read header
	header <- names(read.table(filepath, header = TRUE, sep = sep, nrows = 1, check.names = FALSE))
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
			ordered = list(
				cnvid   = c("cnv_id")
			),
			unordered = list(
				chr   = c("chr", "chrom", "seqnames","chr_child", "chrom_child", "seqnames_child"),
				start = c("start", "begin", "pos_start","start_child", "begin_child", "pos_start_child"),
				end   = c("end", "stop", "pos_end","end_child", "stop_child", "pos_end_child"),
				sample = c("sampleid", "sample_id", "id_sample","sampleid_child"),
				type   = c("type", "svtype", "variant_type","type_child"),
				length = c("size","length","size_child"),
				length_range = c("size_range","length_range"),
				exon_overlap = c("exon_overlap"),
				transcript_overlap = c("transcript_overlap","transcript_bp_overlap"),
				geneid = c("gene_id"),
				overlap_pr = c("overlappr"),
				#biotype = c("biotype"),
				#cnv_inherit = c("inheritance_by_cnv"),
				gene_inherit = c("inheritance_by_gene"),
				loeuf = c("loeuf"),
				overlap_father = c("overlap_cnv_father"),
				overlap_mother = c("overlap_cnv_mother")
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
annotate <- function(cnvs_file, prob_regions_file, output_file, 
										 genome_version = 38, bedtools_path){
	
	gene_annotation_script <- system.file("python",
																				"gene_annotation.py", 
																				package = "MCNV2")
	gene_resource_file <- system.file("resources",
																		"gene_resources.tsv", 
																		package = "MCNV2")
	
	cmd = paste("python3", gene_annotation_script, 
							"--cnv", cnvs_file,
							"--gene_resource", gene_resource_file, 
							"--prob_regions", prob_regions_file,
							"--out", output_file, 
							"--genome_version", genome_version, 
							"--bedtools_path", bedtools_path)
	ret <- system(command = cmd, intern = FALSE)
	
	return(ret)
}


compute_cnv_inheritance <- function(){
	return(NULL)
}

compute_gene_inheritance <- function(){
	return(NULL)
}

#' @export
compute_inheritance <- function(){
	print(system.file("python", "compute_inheritance.py", package = "MCNV2"))
	
	return(0)
}

parse_cnv_size_value <- function(x) {
	if(grepl(x = x, pattern = ">")){
		x <- Inf
	} else {
		x <- gsub("kb", "000", x)
		x <- gsub("Mb", "000000", x)
		x <- as.numeric(x)
	}
	
	return(x)
}

create_dynamic_slider <- function(range_values){
	
	min_val <- range_values[1]
	max_val <- range_values[2]
	diff_val <- max_val - min_val
	
	if (max_val <= 1) {
		# Cas 0–1
		step_val <- 0.1
		seq_vals <- seq(0, 1, by = step_val)
	} else {
		# Cas général : 100 intervalles, arrondis à un chiffre "propre"
		raw_step <- diff_val / 100
		magnitude <- 10 ^ floor(log10(raw_step))
		step_val <- round(raw_step / magnitude) * magnitude
		seq_vals <- seq(
			floor(min_val / step_val) * step_val,
			ceiling(max_val / step_val) * step_val,
			by = step_val
		)
	}
	
	return(seq_vals)
}
