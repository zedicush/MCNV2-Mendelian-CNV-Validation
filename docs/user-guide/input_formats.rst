Input formats
=============

MCNV2 requires two main input files: a **CNV callset** and a **pedigree file**. 

.. important::

   **Column order is critical.** Column names are flexible, but the position of each column must be respected exactly as described below.

CNV callset file
----------------

**Format:** Tab-delimited file

File extension can be ``.txt``, ``.tsv``, ``.bed``, ``.xls``, or any other extension, as long as the content is tab-delimited.

**Required columns (in order):**

1. **CHR** — Chromosome (e.g., chr1, 1, chrX)
2. **START** — CNV start position (integer)
3. **STOP** — CNV end position (integer)
4. **TYPE** — CNV type: ``DEL`` (deletion) or ``DUP`` (duplication)
5. **SAMPLE_ID** — Sample identifier

**Optional columns (after the required 5):**

* **Quality score (highly recommended)** — Caller-specific likelihood score (e.g., PennCNV quality score)
* Number of probes
* Caller concordance percentage (if CNVs are merged from multiple callers)
* Any other metadata

.. important::

   **Quality score (column 6) is highly recommended.**
   
   While optional, the quality score enables generation of **MP versus quality score curves**, which are essential for identifying optimal quality thresholds. Without this column, only size-stratified bar plots will be available.

**Example:**

.. code-block:: text

   chr1    1000000    1050000    DEL    child001    32.5    150
   chr1    2000000    2100000    DUP    child001    45.2    220
   chr2    5000000    5080000    DEL    child002    28.1    95
   
   # Column 6 (32.5, 45.2, 28.1) = Quality scores ← Highly recommended!

.. note::

   * Column names can be anything (e.g., ``chromosome``, ``chr``, ``chrom`` are all acceptable)
   * **Column order cannot change** — position 1 = chromosome, position 2 = start, etc.
   * Optional columns can include quality scores, probe counts, or any other metadata
   * **Caller concordance:** If CNVs from multiple callers are merged, you can include an overlap percentage (e.g., 70% concordance between two callers) as an optional column

Pedigree file
-------------

**Format:** Tab-delimited text file with **exactly 3 columns**

**Required columns (in order):**

1. **SAMPLE_ID** — Offspring/proband sample ID (**must be first column**)
2. **FATHER_ID** — Father sample ID
3. **MOTHER_ID** — Mother sample ID

.. important::

   * **SampleID MUST be in the first column**
   * The order of FatherID and MotherID (columns 2-3) can be swapped if needed
   * The file must have **exactly 3 columns** (no additional columns)
   * All three IDs must form a **complete trio** (no missing parents or offspring)

**Example (correct format):**

.. code-block:: text

   child001    father001    mother001
   child002    father002    mother002
   child003    father003    mother003

**PLINK .fam and KING .kin files:**

If your pedigree file comes from PLINK (``.fam``) or KING (``.kin``), you **must reformat it to 3 columns** before using MCNV2.

PLINK ``.fam`` files typically contain 6 columns:

.. code-block:: text

   FAM001    child001    father001    mother001    1    2
   FAM002    child002    father002    mother002    2    1

**To use with MCNV2, extract only columns 2-4:**

.. code-block:: bash

   # Extract columns 2, 3, 4 from PLINK .fam file
   cut -f2,3,4 pedigree.fam > pedigree_mcnv2.txt

**Result (correct 3-column format):**

.. code-block:: text

   child001    father001    mother001
   child002    father002    mother002

.. note::

   * Column names (headers) are optional
   * Only **complete trios** are analyzed
   * Incomplete trios in the pedigree file are excluded with a warning

Sample ID matching
------------------

**Critical requirement:** The pedigree file must contain **complete trios only**.

**Complete trio definition:**

A trio is complete when the pedigree file contains all three entries:

1. Offspring sample ID
2. Father sample ID  
3. Mother sample ID

.. important::

   **Incomplete trios are excluded:**
   
   * ❌ Father and mother without offspring
   * ❌ Offspring with only one parent
   * ❌ Any trio missing one of the three IDs
   
   **All three IDs must be present in the pedigree file.**

**CNV file presence:**

It is **normal** for a sample to be absent from the CNV file. This simply means the caller did not detect any CNVs for that individual.

**Example scenario:**

.. code-block:: text

   # Pedigree file (complete trio)
   child001    father001    mother001

   # CNV file (only father and child have detected CNVs)
   chr1    1000000    1050000    DEL    child001
   chr2    2000000    2100000    DUP    father001
   # mother001 has NO CNVs detected → this is NORMAL

✅ This trio is **valid** because all three IDs are present in the pedigree file, even though the mother has no detected CNVs.

.. note::

   The absence of CNVs for a sample in the CNV file indicates that the caller found no variants for that individual. This does not invalidate the trio.

Reference annotations
---------------------

MCNV2 uses the following reference files (provided with the package):

**Problematic regions** (UCSC-based, BED format):

* Segmental duplications
* Centromeres
* Telomeres
* HLA region

**Gene annotations** (tab-delimited or BED):

* Gene coordinates (from GTF/GFF)
* LOEUF scores (gnomAD v4) — loss-of-function constraint metric

**Gene exclusion list** (optional, tab-delimited):

* Ensembl gene IDs to exclude from Mendelian Precision calculation

.. note::

   **Genome build:** hg38 is supported in this release. hg19/GRCh37 support is planned.

File validation
---------------

Before running MCNV2, verify:

1. ✅ CNV file has correct column order (chr, start, stop, type, sample_id, optional...)
2. ✅ Pedigree file has **exactly 3 columns** (offspring, father, mother)
3. ✅ Pedigree file contains **only complete trios** (all 3 IDs present for each trio)
4. ✅ CNV file is tab-delimited (not comma-separated)
5. ✅ Chromosome names are consistent (e.g., always ``chr1`` or always ``1``)
6. ✅ CNV types are ``DEL`` or ``DUP`` (case-insensitive)

.. tip::

   Use the Shiny interface to upload files — MCNV2 will automatically validate formats and report any issues.
