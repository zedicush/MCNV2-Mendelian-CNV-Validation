Overview
========
MCNV2 (Mendelian CNV Validation) is a framework designed to assess the quality of copy number variations (CNVs) using **family-based Mendelian inheritance** in parent–offspring trios.

MCNV2 is designed for: - **Researchers** optimizing CNV calling pipelines (WGS/WES/arrays) - **Bioinformatics pipelines** that require reproducible QC summaries - **Clinical and research labs** comparing callers or thresholds

--------------

Key idea
--------

| True de novo CNVs are rare (1.92%), so most CNVs observed only in offspring are enriched for technical false positives.
| MCNV2 leverages this signal to compute **Mendelian Precision (MP)** without requiring an external benchmark callset.

--------------

What MCNV2 provides
-------------------

-  **CLI** for automated analyses (pipeline integration, large cohorts)
-  **Shiny app** for interactive exploration (filters + plots + exports)
-  **Annotation** modules:

   -  gene overlap
   -  gene constraint (LOEUF; gnomAD)
   -  problematic regions (UCSC tracks: segdups, centromeres, telomeres, HLA)

-  **Inheritance modes**

   -  CNV-level (reciprocal overlap threshold)
   -  gene-level (shared affected genes, robust to breakpoint noise)

-  **Exports**

   -  MP summary tables
   -  candidate de novo / Mendelian errors lists

--------------

Recommended starting workflow
-----------------------------

1. Install MCNV2 and dependencies
2. Run preprocessing/annotation
3. Compute MP with baseline filtering
4. Compare filtering strategies (score / concordance / LOEUF)
5. Export high-confidence CNVs and de novo candidates
