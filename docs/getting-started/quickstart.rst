Quickstart
==========

This page provides a high-level overview of the MCNV2 workflow, from data input to interpretation of Mendelian Precision results.

MCNV2 can be used either through an interactive Shiny application or via a command-line interface (CLI) for batch and pipeline-based analyses.

--------------

Basic workflow
--------------

Step 1 — Provide input data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Upload the following files:

-  **CNV file** (tab-delimited): chromosome, start, end, CNV type (DEL/DUP), sample ID
-  **Pedigree file**: parent–offspring relationships (PLINK ``.fam`` or KING ``.kin``)

Only complete trios (child + both parents) are retained for analysis.

--------------

Step 2 — Choose inheritance mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select how Mendelian inheritance is evaluated:

-  | **CNV-level (default)**
   | Inheritance is defined using reciprocal genomic overlap between child and parental CNVs.

-  | **Gene-level**
   | Inheritance is defined at the gene level, allowing robustness to breakpoint imprecision.

--------------

Step 3 — Set overlap parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Reciprocal overlap threshold**
   Default: **50%**
   Adjustable from 1% to 100%, depending on desired stringency.

--------------

Step 4 — Run Mendelian Precision analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Click **“Go to Mendelian Precision analysis”** to compute Mendelian Precision (MP) across:

-  CNV type (deletions vs duplications)
-  CNV size ranges
-  Quality metrics (if available)
-  Optional filters

--------------

Running MCNV2 via CLI
---------------------

MCNV2 can also be run in batch mode using the command-line interface (CLI), making it suitable for large cohorts and reproducible pipelines.

!!! note CLI commands depend on your local installation and wrapper scripts. See the full CLI tutorial for recommended usage patterns and examples.

👉 **Go to:** [CLI tutorial -> (../tutorials/cli_tutorial.md)

--------------

What to look at first
---------------------

When exploring results, we recommend the following progression:

1. Baseline Mendelian Precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  MP computed with minimal filtering
-  Provides a global view of call quality

2. MP versus quality score
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Size-stratified curves
-  Identify score thresholds where MP plateaus

3. Impact of key filters
~~~~~~~~~~~~~~~~~~~~~~~~

Evaluate how MP changes when applying:

-  **Caller concordance** (≥2 algorithms)
-  **LOEUF filtering** (exclude highly constrained genes)
-  **Problematic region exclusion** (segmental duplications, centromeres, HLA, etc.)

--------------

Outputs
-------

MCNV2 produces the following outputs:

-  **Mendelian Precision summary tables**

   -  DEL vs DUP
   -  CNV size bins
   -  Filtering configurations

-  **Publication-ready figures**

   -  Exportable from the Shiny interface

-  **Lists of Mendelian errors**

   -  Candidate *de novo* CNVs
   -  Residual technical artifacts for downstream validation

All tables and figures can be downloaded for further analysis or reporting.
