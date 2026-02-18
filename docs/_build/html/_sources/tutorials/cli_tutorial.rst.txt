CLI tutorial
============

MCNV2 provides R functions for **reproducible batch processing** and integration into automated pipelines.

.. note::

   These functions are R wrappers that call Python scripts internally. You must have 
   the MCNV2 R package installed and Python environment configured 
   (see :doc:`../getting-started/installation`).

Overview
--------

MCNV2 provides three main functions for command-line workflows:

1. **annotate()** — Annotate CNVs with genes, LOEUF scores, and problematic regions
2. **compute_inheritance()** — Calculate inheritance status (Transmitted_CNV and Transmitted_gene)
3. **compute_mp()** — Calculate Mendelian Precision with filtering options

Typical workflow
----------------

.. code-block:: R

   library(MCNV2)
   
   # 1. Annotate CNVs
   annotate(
     cnvs_file = "data/cnvs.tsv",
     prob_regions_file = "data/problematic_regions.bed",
     output_file = "results/cnvs_annotated.tsv",
     genome_version = 38,
     bedtools_path = "/usr/local/bin/bedtools"
   )
   
   # 2. Compute inheritance status
   compute_inheritance(
     cnvs_file = "results/cnvs_annotated.tsv",
     pedigree_file = "data/pedigree.tsv",
     output_file = "results/cnvs_inheritance.tsv",
     overlap = 0.5
   )
   
   # 3. Compute Mendelian Precision
   # Note: MP is stratified by CNV type (DEL vs DUP) by default
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_summary.tsv",
     transmission_type = "cnv",
     min_size = 30000,
     max_prob_regions = 0.5
   )

Function reference
------------------

annotate()
~~~~~~~~~~

Annotate CNVs with gene information, LOEUF scores, and problematic region overlap.

**Usage:**

.. code-block:: R

   annotate(
     cnvs_file,
     prob_regions_file,
     output_file,
     genome_version = 38,
     bedtools_path
   )

**Parameters:**

* **cnvs_file** (character) — Path to CNV file (tab-delimited)

  Required columns: CHR, START, STOP, TYPE, SAMPLE_ID

* **prob_regions_file** (character) — Path to problematic regions BED file

  Default file provided in package if not specified

* **output_file** (character) — Path for output annotated file

* **genome_version** (numeric) — Genome build: 38 (GRCh38/hg38) or 37 (GRCh37/hg19)

  Default: 38

* **bedtools_path** (character) — Path to bedtools executable

  Example: "/usr/local/bin/bedtools"

**Returns:**

* 0 if successful
* 1 if output file was not created (error)

**Output columns:**

All input columns plus:

* GeneName — HGNC gene symbol
* GeneID — Ensembl gene ID
* Transcript — Ensembl transcript ID
* LOEUF — Loss-of-function constraint score (gnomAD v4)
* problematic_region_overlap — Percentage overlap with problematic regions

**Example:**

.. code-block:: R

   library(MCNV2)
   
   # Annotate CNVs
   status <- annotate(
     cnvs_file = "data/cnvs.tsv",
     prob_regions_file = system.file("resources", "problematic_regions.bed", 
                                     package = "MCNV2"),
     output_file = "results/cnvs_annotated.tsv",
     genome_version = 38,
     bedtools_path = "/usr/local/bin/bedtools"
   )
   
   if (status == 0) {
     message("Annotation completed successfully")
   } else {
     stop("Annotation failed")
   }

compute_inheritance()
~~~~~~~~~~~~~~~~~~~~~

Calculate inheritance status for each CNV based on parental data.

**Usage:**

.. code-block:: R

   compute_inheritance(
     cnvs_file,
     pedigree_file,
     output_file,
     overlap = 0.5
   )

**Parameters:**

* **cnvs_file** (character) — Path to annotated CNV file (output from annotate())

* **pedigree_file** (character) — Path to pedigree file

  Three columns: SAMPLE_ID, FATHER_ID, MOTHER_ID (tab-delimited, no header)

* **output_file** (character) — Path for output file with inheritance status

* **overlap** (numeric) — Minimum overlap for a CNV to be considered as inherited

  Range: 0.01 to 1.0
  
  Default: 0.5 (50%)

**Returns:**

* 0 if successful
* 1 if output file was not created (error)

**Output columns:**

All input columns plus:

* **Transmitted_CNV** — True/False (coordinate-based inheritance)
* **Transmitted_gene** — True/False/intergenic (gene-based inheritance)

**Example:**

.. code-block:: R

   library(MCNV2)
   
   # Compute inheritance with 50% overlap threshold
   status <- compute_inheritance(
     cnvs_file = "results/cnvs_annotated.tsv",
     pedigree_file = "data/pedigree.tsv",
     output_file = "results/cnvs_inheritance.tsv",
     overlap = 0.5
   )
   
   if (status == 0) {
     message("Inheritance calculation completed")
     
     # Read results
     cnvs <- read.table("results/cnvs_inheritance.tsv", 
                        header = TRUE, sep = "\t", quote = "")
     
     # Count inherited vs non-inherited
     table(cnvs$Transmitted_CNV)
   }

compute_mp()
~~~~~~~~~~~~

Calculate Mendelian Precision with optional filtering.

.. important::

   **MP is always stratified by CNV type (DEL vs DUP) by default.** 
   
   This is the recommended approach as deletions and duplications have different 
   quality profiles. Set ``stratify_by_type = FALSE`` only if you specifically need 
   a single global MP value (not recommended).

.. important::

   This function requires the output file from ``compute_inheritance()`` as input.

**Usage:**

.. code-block:: R

   compute_mp(
     inheritance_file,
     output_file,
     transmission_type = "cnv",
     min_size = NULL,
     max_size = NULL,
     min_score = NULL,
     max_prob_regions = NULL,
     min_loeuf = NULL,
     stratify_by_size = FALSE,
     stratify_by_type = TRUE
   )

**Parameters:**

* **inheritance_file** (character) — Path to inheritance file

  **Must be the output from compute_inheritance()**
  
  Required columns: Transmitted_CNV, Transmitted_gene

* **output_file** (character) — Path for MP summary output

* **transmission_type** (character) — Transmission matching type

  * "cnv" — Use Transmitted_CNV (coordinate-based)
  * "gene" — Use Transmitted_gene (gene-based)
  
  Default: "cnv"

* **min_size** (numeric) — Minimum CNV size in bp (optional)

  Example: 30000 (30 kb)

* **max_size** (numeric) — Maximum CNV size in bp (optional)

* **min_score** (numeric) — Minimum quality score (optional)

* **max_prob_regions** (numeric) — Maximum problematic regions overlap (0-1, optional)

  Example: 0.5 (exclude CNVs with >50% overlap)

* **min_loeuf** (numeric) — Minimum LOEUF threshold (optional)

  Example: 0.6 (exclude constrained genes)

* **stratify_by_size** (logical) — Stratify MP by size ranges

  Default: FALSE

* **stratify_by_type** (logical) — Stratify MP by CNV type (DEL/DUP)

  Default: TRUE (recommended)
  
  **Keep TRUE** unless you specifically need a single global MP value. 
  Deletions and duplications have different quality profiles and should be 
  evaluated separately.

**Returns:**

* 0 if successful
* 1 if output file was not created (error)

**Output:**

Tab-delimited file with MP statistics:

* CNV_type (DEL/DUP, or All if stratify_by_type=FALSE)
* Size_range (All, or specific ranges if stratify_by_size=TRUE)
* Total_CNVs
* Inherited_CNVs
* Non_inherited_CNVs
* MP (Mendelian Precision %)

**Example output (default: stratified by type):**

.. code-block:: text

   CNV_type  Size_range  Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   DEL       All         5000        4250            750                  85.0
   DUP       All         3000        2400            600                  80.0

**Example output (stratified by type AND size):**

.. code-block:: text

   CNV_type  Size_range   Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   DEL       1-30kb       800         600             200                  75.0
   DEL       30-50kb      600         540             60                   90.0
   DEL       50-100kb     900         855             45                   95.0
   DUP       1-30kb       600         420             180                  70.0
   DUP       30-50kb      500         425             75                   85.0
   DUP       50-100kb     700         665             35                   95.0

**Example 1: Basic MP calculation (stratified by type)**

.. code-block:: R

   library(MCNV2)
   
   # Calculate MP (CNV-level, stratified by DEL/DUP, no filters)
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_summary.tsv",
     transmission_type = "cnv"
   )
   
   # Read MP results
   mp <- read.table("results/mp_summary.tsv", header = TRUE, sep = "\t")
   print(mp)
   #   CNV_type  Size_range  Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   #   DEL       All         5000        4250            750                  85.0
   #   DUP       All         3000        2400            600                  80.0

**Example 2: MP with filtering**

.. code-block:: R

   # Calculate MP with size and quality filters
   # Still stratified by DEL/DUP (default)
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_filtered.tsv",
     transmission_type = "cnv",
     min_size = 30000,           # ≥30 kb
     max_prob_regions = 0.5,     # ≤50% prob regions overlap
     min_score = 100             # Score ≥100
   )
   
   # Output shows DEL and DUP separately after filtering
   #   CNV_type  Size_range  Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   #   DEL       All         3000        2775            225                  92.5
   #   DUP       All         1800        1620            180                  90.0

**Example 3: Technical MP (excluding constrained genes)**

.. code-block:: R

   # Calculate technical MP (excluding LOEUF < 0.6)
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_technical.tsv",
     transmission_type = "gene",
     min_loeuf = 0.6
   )
   
   # Output:
   #   CNV_type  Size_range  Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   #   DEL       All         4500        4185            315                  93.0
   #   DUP       All         2700        2511            189                  93.0

**Example 4: MP stratified by size**

.. code-block:: R

   # Calculate MP for each size range
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_by_size.tsv",
     transmission_type = "cnv",
     stratify_by_size = TRUE,
     stratify_by_type = TRUE
   )
   
   # Output shows DEL and DUP for each size range
   #   CNV_type  Size_range   Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   #   DEL       1-30kb       800         600             200                  75.0
   #   DEL       30-50kb      600         540             60                   90.0
   #   DEL       50-100kb     900         855             45                   95.0
   #   DEL       100-200kb    700         686             14                   98.0
   #   DUP       1-30kb       600         420             180                  70.0
   #   DUP       30-50kb      500         425             75                   85.0
   #   DUP       50-100kb     700         665             35                   95.0
   #   DUP       100-200kb    500         490             10                   98.0

**Example 5: Global MP (not recommended)**

.. code-block:: R

   # Calculate single global MP value (not stratified by type)
   # NOT RECOMMENDED: loses information about DEL vs DUP differences
   compute_mp(
     inheritance_file = "results/cnvs_inheritance.tsv",
     output_file = "results/mp_global.tsv",
     transmission_type = "cnv",
     stratify_by_type = FALSE
   )
   
   # Output (single row):
   #   CNV_type  Size_range  Total_CNVs  Inherited_CNVs  Non_inherited_CNVs   MP
   #   All       All         8000        6650            1350                 83.1

Batch processing example
------------------------

Process multiple datasets:

.. code-block:: R

   library(MCNV2)
   
   # List of datasets
   datasets <- c("cohort1", "cohort2", "cohort3")
   
   for (dataset in datasets) {
     message(paste("Processing", dataset))
     
     # 1. Annotate
     annotate(
       cnvs_file = paste0("data/", dataset, "_cnvs.tsv"),
       prob_regions_file = system.file("resources", "problematic_regions.bed", 
                                       package = "MCNV2"),
       output_file = paste0("results/", dataset, "_annotated.tsv"),
       genome_version = 38,
       bedtools_path = "/usr/local/bin/bedtools"
     )
     
     # 2. Inheritance
     compute_inheritance(
       cnvs_file = paste0("results/", dataset, "_annotated.tsv"),
       pedigree_file = paste0("data/", dataset, "_pedigree.tsv"),
       output_file = paste0("results/", dataset, "_inheritance.tsv"),
       overlap = 0.5
     )
     
     # 3. MP (multiple strategies)
     # All MP calculations are stratified by type (DEL/DUP) by default
     
     # 3a. Overall MP
     compute_mp(
       inheritance_file = paste0("results/", dataset, "_inheritance.tsv"),
       output_file = paste0("results/", dataset, "_mp_all.tsv"),
       transmission_type = "cnv"
     )
     
     # 3b. Filtered MP
     compute_mp(
       inheritance_file = paste0("results/", dataset, "_inheritance.tsv"),
       output_file = paste0("results/", dataset, "_mp_filtered.tsv"),
       transmission_type = "cnv",
       min_size = 30000,
       max_prob_regions = 0.5
     )
     
     # 3c. Technical MP
     compute_mp(
       inheritance_file = paste0("results/", dataset, "_inheritance.tsv"),
       output_file = paste0("results/", dataset, "_mp_technical.tsv"),
       transmission_type = "gene",
       min_loeuf = 0.6
     )
   }
   
   # Combine results (each file has DEL and DUP rows)
   all_mp <- do.call(rbind, lapply(datasets, function(d) {
     mp <- read.table(paste0("results/", d, "_mp_all.tsv"), 
                      header = TRUE, sep = "\t")
     mp$Dataset <- d
     return(mp)
   }))
   
   # Result has DEL and DUP rows for each dataset
   write.table(all_mp, "results/combined_mp.tsv", 
               sep = "\t", row.names = FALSE, quote = FALSE)

Pipeline integration
--------------------

Nextflow workflow
~~~~~~~~~~~~~~~~~

.. code-block:: groovy

   // main.nf
   process annotate {
       input:
       path cnv
       path prob_regions
       
       output:
       path "${cnv.baseName}_annotated.tsv"
       
       script:
       """
       Rscript -e '
       MCNV2::annotate(
           cnvs_file = "${cnv}",
           prob_regions_file = "${prob_regions}",
           output_file = "${cnv.baseName}_annotated.tsv",
           genome_version = 38,
           bedtools_path = "/usr/local/bin/bedtools"
       )'
       """
   }
   
   process inheritance {
       input:
       path annotated
       path pedigree
       
       output:
       path "${annotated.baseName}_inheritance.tsv"
       
       script:
       """
       Rscript -e '
       MCNV2::compute_inheritance(
           cnvs_file = "${annotated}",
           pedigree_file = "${pedigree}",
           output_file = "${annotated.baseName}_inheritance.tsv",
           overlap = 0.5
       )'
       """
   }
   
   process mp {
       input:
       path inheritance
       
       output:
       path "${inheritance.baseName}_mp.tsv"
       
       script:
       """
       Rscript -e '
       MCNV2::compute_mp(
           inheritance_file = "${inheritance}",
           output_file = "${inheritance.baseName}_mp.tsv",
           transmission_type = "cnv",
           min_size = 30000
       )'
       """
   }

Error handling
--------------

Check return codes:

.. code-block:: R

   library(MCNV2)
   
   # Annotate with error checking
   status <- annotate(
     cnvs_file = "data/cnvs.tsv",
     prob_regions_file = "data/problematic_regions.bed",
     output_file = "results/annotated.tsv",
     genome_version = 38,
     bedtools_path = "/usr/local/bin/bedtools"
   )
   
   if (status != 0) {
     stop("Annotation failed. Check input files and bedtools path.")
   }
   
   # Verify output exists and is not empty
   if (!file.exists("results/annotated.tsv")) {
     stop("Output file was not created")
   }
   
   if (file.size("results/annotated.tsv") == 0) {
     stop("Output file is empty")
   }
   
   message("Annotation completed successfully")

Performance considerations
--------------------------

**Memory requirements:**

* Annotation: ~4 GB for 100,000 CNVs
* Inheritance: ~8 GB for 100,000 CNVs in 1,000 trios
* MP calculation: ~2 GB

**Runtime (approximate):**

* Annotation: ~5 min for 100,000 CNVs
* Inheritance: ~10 min for 100,000 CNVs in 1,000 trios
* MP calculation: ~1 min

**Recommendations:**

* For large datasets (>500,000 CNVs), consider splitting by chromosome
* Use parallel processing for batch analysis of multiple cohorts
* Ensure sufficient disk space for intermediate files

See also
--------

* :doc:`../user-guide/preprocessing` — Preprocessing steps explained
* :doc:`../user-guide/inheritance` — Inheritance calculation details
* :doc:`../user-guide/mendelian_precision` — MP calculation and interpretation
* :doc:`../user-guide/filtering` — Filtering strategies
