Preprocessing
=============

The **Preprocessing** tab is the first step in the MCNV2 workflow. It performs CNV annotation and inheritance status calculation.

Overview
--------

This tab allows you to:

1. Upload input files (CNV calls, pedigree, problematic regions)
2. Set inheritance parameters (overlap threshold, genome build)
3. Annotate CNVs with genes, LOEUF scores, and problematic regions
4. Calculate inheritance status (transmitted vs non-transmitted)
5. View annotated results

Input files
-----------

**CNV file** (mandatory)

Tab-delimited file with columns: CHR, START, STOP, TYPE (DEL/DUP), SAMPLE_ID

See :doc:`input_formats` for detailed specifications.

**Pedigree file** (mandatory)

Three-column file: SAMPLE_ID, FATHER_ID, MOTHER_ID

Only complete trios are analyzed.

**Problematic regions file** (optional, BED format)

BED file with problematic genomic regions:

* Segmental duplications
* Centromeres
* Telomeres
* HLA region

.. note::

   A default file is provided. You can replace it with your own BED file if needed.

Parameters for inheritance calculation
--------------------------------------

**Inheritance threshold (child CNV proportion)**

Default: **0.5** (50%)

This parameter defines the minimum reciprocal overlap required between a child CNV and a parental CNV to consider the child CNV as inherited.

* **Value range:** 0.01 to 1.0 (1% to 100%)
* **Interpretation:** If ≥X% of a child CNV overlaps with a parental CNV (father or mother), the CNV is classified as **inherited (True)**
* **Example:** With threshold=0.5, a child CNV is inherited if at least 50% of it overlaps with a CNV in at least one parent

**Genome build**

* **Default:** GRCh38/hg38
* **Planned:** GRCh37/hg19 (not yet implemented)

Workflow
--------

**Step 1: Upload files**

Upload your CNV file, pedigree file, and optionally a custom problematic regions file.

**Step 2: Set parameters**

* **Inheritance threshold:** Adjust the reciprocal overlap percentage (default 50%)
* **Genome build:** Select GRCh38/hg38

**Step 3: Submit**

Click **Submit** to start annotation and inheritance calculation.

Annotation process
------------------

MCNV2 annotates each CNV with:

**Gene annotation**

Each CNV is intersected with gene coordinates (Gencode v45):

* **GeneName** — HGNC gene symbol
* **GeneID** — Ensembl gene ID
* **Transcript** — Ensembl transcript ID

If a CNV overlaps multiple genes, one row is created per gene (CNV-gene pairs).

**LOEUF scores** (gnomAD v4)

LOEUF (Loss-of-function Observed/Expected Upper bound Fraction) quantifies gene constraint:

* **Low LOEUF** (≤0.6): Highly constrained genes (intolerant to loss-of-function)
* **High LOEUF** (>0.6): Less constrained genes

LOEUF is used for:

* Stratifying Mendelian Precision by gene constraint
* Optional filtering (exclude constrained genes to focus on technical precision)

**Problematic regions overlap**

Percentage of CNV overlapping with problematic regions:

* Segmental duplications
* Centromeres
* Telomeres
* HLA region

.. note::

   This percentage is used in the filtering step to exclude CNVs with high overlap (e.g., >50%).

Annotated CNV table
-------------------

After clicking **Submit**, the first table displays annotated CNVs:

**Columns:**

* Original CNV file columns (CHR, START, STOP, TYPE, SAMPLE_ID, quality scores, etc.)
* **GeneName** — Overlapping gene name
* **GeneID** — Ensembl gene ID
* **Transcript** — Ensembl transcript ID
* **LOEUF** — Constraint score
* **problematic_region_overlap** — Percentage overlap with problematic regions

**Table navigation:**

* Show 10, 50, 100, or all entries
* Scroll horizontally to view all columns

**File path:**

The file path where the annotated table is saved is displayed below the table.

Inheritance status calculation
-------------------------------

Click **Proceed to inheritance status** to calculate transmission for each CNV.

MCNV2 uses two complementary approaches to determine inheritance:

1. **Transmitted_CNV** (coordinate-based) — Based on genomic coordinate overlap
2. **Transmitted_gene** (gene-based) — Based on shared affected genes

.. seealso::

   See :doc:`inheritance` for a comprehensive explanation of the two inheritance matching approaches.

**Column values in the inheritance table:**

* **Transmitted_CNV:** True / False
* **Transmitted_gene:** True / False / intergenic

Inheritance status table
------------------------

The second table displays inheritance results:

**Columns:**

* All columns from the annotated CNV table
* **Transmitted_CNV** — True/False (coordinate-based inheritance)
* **Transmitted_gene** — True/False/intergenic (gene-based inheritance)

**Interpretation:**

* **True (both columns):** CNV is inherited from at least one parent
* **False (Transmitted_CNV):** Candidate *de novo* CNV (no parental overlap)
* **False (Transmitted_gene):** Gene not affected in parents
* **"intergenic" (Transmitted_gene):** CNV does not overlap any gene

**File path:**

The file path where the inheritance table is saved is displayed below the table.

Next step: Mendelian Precision analysis
----------------------------------------

Once inheritance status is calculated, click **Go to Mendelian Precision analysis** to:

* Compute Mendelian Precision across size ranges and quality thresholds
* Apply filters (quality scores, problematic regions, LOEUF, caller concordance)
* Generate publication-ready plots

See :doc:`mendelian_precision` for details on the MP analysis workflow.

Tips
----

**Understanding the two inheritance approaches**

Both **Transmitted_CNV** (coordinate-based) and **Transmitted_gene** (gene-based) have different strengths:

* **Transmitted_CNV** requires an overlap threshold but works for all CNVs
* **Transmitted_gene** is more robust to breakpoint and CNV fragmentation but only works for genic CNVs

**Recommendation:** Use both approaches to get a comprehensive view of CNV inheritance.

See :doc:`inheritance` for detailed comparison, advantages, limitations, and use cases.

**Filtering intergenic CNVs**

Intergenic CNVs can only be evaluated using coordinate-based matching (Transmitted_CNV).

**Table too large to display**

Use pagination (show 10/50/100 entries).

Troubleshooting
---------------

**Error: "No complete trios found"**

Check that:

* Pedigree file has correct format (SAMPLE_ID, FATHER_ID, MOTHER_ID)
* All three IDs (child, father, mother) are present in the pedigree file
* Sample IDs are consistent between CNV and pedigree files

**Warning: "Some samples in pedigree not found in CNV file"**

This is normal. If a sample has no detected CNVs, it won't appear in the CNV file. The trio is still valid.

**Table too large to display**

Use pagination (show 10/50/100 entries) or download the full table for offline analysis.
