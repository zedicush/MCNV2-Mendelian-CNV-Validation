Mendelian Precision
===================

Mendelian Precision (MP) is the core quality metric in MCNV2, quantifying the proportion of CNV calls that follow Mendelian inheritance patterns in parent–offspring trios.

.. note:: 
  
   Mendelian Precision is defined as:

   .. math::

      \text{MP} = \frac{I}{N} = 1 - \frac{E}{N}

   where:

   * **N** = Total number of CNVs detected in offspring
   * **I** = Number of inherited CNVs (found in at least one parent)
   * **E** = Number of non-inherited CNVs (not found in parents)

**Rationale:**

Given the low expected rate of genuine *de novo* CNVs (typically 1.92% per individual), non-inherited CNVs are predominantly false positives. MP therefore provides a biologically grounded estimate of CNV call precision without requiring an external reference callset.

Interpretation
--------------

**High MP (≥80%):**

* Most CNVs are inherited → High call quality
* Few false positives
* Callset suitable for downstream analyses

**Moderate MP (50-80%):**

* Mixed quality
* May benefit from additional filtering
* Investigate size-specific or quality-score-specific patterns

**Low MP (<50%):**

* High false positive rate
* Requires aggressive filtering or caller parameter tuning
* Consider alternative CNV calling methods

Transmission types
------------------

MCNV2 calculates MP using two complementary approaches:

1. **CNV-level matching** (coordinate-based) — Based on genomic coordinate overlap
2. **Gene-level matching** (gene-based) — Based on shared affected genes

.. seealso::

   See :doc:`inheritance` for detailed algorithms, advantages, and limitations 
   of each approach.

**MP calculation:**

MP can be computed separately using each approach:

.. math::

   \text{MP}_{\text{CNV}} = \frac{\text{True (Transmitted_CNV)}}{\text{Total CNVs}}

.. math::  

   \text{MP}_{\text{gene}} = \frac{\text{True (Transmitted_gene)}}{\text{Total genic CNVs}}

**Why use both:**

* CNV-level works for all CNVs (including intergenic)
* Gene-level is robust to breakpoint variability (genic CNVs only)
* The two approaches provide complementary quality assessments

**Recommendation:** 

Compute MP using both approaches. Differences between MP_CNV and MP_gene 
reflect the interplay of breakpoint precision and CNV fragmentation in 
your specific dataset.

MP stratification
-----------------

To identify quality issues and optimization strategies, MP should be computed across multiple dimensions:

By CNV size
~~~~~~~~~~~

**Size ranges:**

* 1-30 kb
* 30-50 kb
* 50-100 kb
* 100-200 kb
* 200-500 kb
* 500 kb-1 Mb
* >1 Mb

**Typical pattern:**

* Small CNVs (<30 kb): Low MP (high false positive rate)
* Medium CNVs (50-200 kb): Moderate to high MP
* Large CNVs (>500 kb): High MP (more reliable)

**Interpretation:**

If MP increases dramatically with size, consider applying a minimum size filter.

By CNV type
~~~~~~~~~~~

**DEL vs DUP:**

MP is calculated separately for deletions and duplications, as they often have different precision profiles.

By quality score
~~~~~~~~~~~~~~~~

**Quality metrics:**

MP can be stratified by various quality scores:

* **Score** — Caller-specific quality score
* **SNP** — Number of supporting SNPs (array data)
* **% Overlap** — Reciprocal overlap percentage between CNVs detected by multiple algorithms

**Typical pattern:**

MP increases with quality score, often plateauing at a certain threshold. This plateau identifies the optimal filtering threshold.

**Example:**

* Score ≥10: MP = 65%
* Score ≥30: MP  780%
* Score ≥50: MP= 90%
* Score ≥200: M = 9 2% (plateau)

**Optimal threshold:** 50 (additional filtering beyond this point provides no MP improvement)

Filtering strategies
--------------------

To maximize MP while retaining biologically relevant CNVs, MCNV2 supports multiple filtering approaches:

CNV-level filters
~~~~~~~~~~~~~~~~~

**Size filters:**

* Minimum size (bp)
* Maximum size (bp)
* Target specific size ranges

**Quality score filters:**

* Minimum quality score
* Minimum number of supporting probes/reads
* Minimum caller concordance

.. important::

   **Highly recommended:** Apply problematic region filters to all datasets. 
   
   CNVs in these regions are enriched for false positives due to:
   
   * Read mismapping (segmental duplications)
   * Poor mappability (centromeres, telomeres)
   * Genuine polymorphism (HLA region)
   
   Recommended threshold: Exclude CNVs with >50% overlap.

Gene-level filters
~~~~~~~~~~~~~~~~~~

**Exclusion lists:**

* Upload a list of genes to exclude (e.g., immunoglobulin genes, olfactory receptors)
* Useful for removing genes prone to technical artifacts

**Gene constraint filters (LOEUF):**

* Exclude CNVs affecting highly constrained genes (LOEUF < threshold)
* Rationale: CNVs in constrained genes may be genuine *de novo* events, reducing technical MP

.. important::

   **LOEUF filter for MP calculation only**
   
   Excluding CNVs in constrained genes (low LOEUF) helps distinguish:
   
   * **Technical MP** — Precision excluding likely *de novo* events
   * **Overall MP** — Precision including all non-inherited CNVs
   
   CNVs affecting constrained genes are enriched for genuine *de novo* events, 
   which reduce MP but are biologically real. Excluding them from MP calculation 
   provides a cleaner assessment of **technical false positive rate**.
   
   **Critical:** These CNVs must be **retained** in your final dataset for 
   downstream analyses (disease association, burden tests) as they may represent 
   pathogenic variants.

See :doc:`filtering` for detailed filtering strategies and optimization approaches.

Optimal threshold identification
---------------------------------

**Strategy:**

1. Plot MP versus quality score threshold (line plot)
2. Identify where MP plateaus
3. Use the lowest threshold at which MP reaches plateau

**Example:**

.. code-block:: text

   For deletions 50-100kb:
   - Score ≥30  → MP = 75% (n = 500)
   - Score ≥50 → MP = 85% (n = 300)
   - Score ≥70 → MP = 92% (n = 150)  ← Plateau starts
   - Score ≥70 → MP = 92% (n = 140)  ← No further MP gain
   
   Optimal threshold: 150 (plateau reached, retains 150 CNVs)


Technical vs biological non-inheritance
---------------------------------------

**Challenge:**

Non-inherited CNVs include both:

1. **Technical false positives** (reduce MP, should be filtered)
2. **Genuine *de novo* CNVs** (reduce MP, should be retained)

**Solution:**

Use gene constraint (LOEUF) to distinguish these two categories:

* CNVs affecting constrained genes (low LOEUF) → Enriched for *de novo* events
* CNVs affecting unconstrained genes → Enriched for false positives

**Workflow:**

1. Calculate **MP (all CNVs)** → Includes both technical and biological non-inheritance
2. Calculate **MP (excluding LOEUF < 0.6)** → Focuses on technical precision only
3. Compare the two values:

   * **Large difference (>10%):** High *de novo* rate in constrained genes
   * **Small difference (<5%):** Most non-inherited CNVs are false positives

.. important::

   **LOEUF exclusion is for assessment only**
   
**Example:**

.. code-block:: text

   1000 CNVs total:
   - MP (all CNVs) = 85% → 150 non-inherited
   - MP (LOEUF ≥ 0.6) = 92% → 80 non-inherited (970 CNVs after excluding 30)
   
   Interpretation:
   - ~30 CNVs (3%) likely genuine *de novo* in constrained genes
   - ~120 CNVs (12%) technical false positives
   - Keep all 1000 CNVs for downstream analysis
   - Apply quality filters to remove the 120 false positives


Global vs stratified MP
------------------------

**Global MP:**

* Single value for all CNVs (or all DEL, all DUP)
* Simple summary metric
* May mask quality issues in specific size ranges or quality bins

**Stratified MP:**

* MP computed for each size range and/or quality threshold
* Reveals patterns (e.g., low MP for small CNVs)
* Enables targeted filtering strategies

**Best practice:** Always examine stratified MP before applying filters.

.. admonition:: Use cases for Mendelian Precision

  **1. CNV caller evaluation**

  Compare MP across different callers or parameter settings to identify the best configuration.

  **2. Quality control**

  Assess CNV call quality in a new dataset.

  **3. Filter optimization**

  Systematically evaluate the impact of different filters on MP to maximize precision while retaining CNVs.

  **4. Method development**

  Use MP as an optimization criterion when developing new CNV calling methods.

  **5. Publication**

  Report MP alongside standard metrics (sensitivity, specificity) to provide a biologically grounded quality assessment.

Limitations
-----------

**Assumptions:**

* Parents are biological parents (not step-parents or adoption)
* Trio relationships are correctly specified in pedigree file
* *De novo* CNV rate is low (~1.92%)

**Caveats:**

* MP does not distinguish between technical false positives and genuine *de novo* events
* MP is not a measure of sensitivity (cannot detect false negatives)

See also
--------

* :doc:`preprocessing` — How inheritance status is calculated
* :doc:`inheritance` — CNV-level vs gene-level matching in detail
* :doc:`filtering` — Filtering strategies to optimize MP
* :doc:`outputs` — Visualizations and tables
