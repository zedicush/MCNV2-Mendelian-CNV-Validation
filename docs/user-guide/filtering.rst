Filtering strategies
====================

Filtering is essential to maximize Mendelian Precision while retaining biologically relevant CNVs. MCNV2 provides flexible filtering options across multiple dimensions.

Overview
--------

**Goal of filtering:**

* Increase Mendelian Precision by removing false positives
* Balance precision and sensitivity

.. admonition:: **Two-stage approach:**

   1. **MP Exploration** — Apply broad filters to assess baseline quality
   2. **Fine-tuning** — Refine filters systematically to identify optimal thresholds

**Key principles:**

* Calculate MP separately for deletions and duplications
* Stratify by CNV size (small vs large CNVs have different quality profiles)
* Consider quality metrics
* Consider gene constraint when evaluating MP

Filter categories
-----------------

Size-based filters
~~~~~~~~~~~~~~~~~~

**CNV size (bp):**

* **Minimum size:** Exclude very small CNVs (e.g., <10 kb)
* **Maximum size:** Optionally exclude very large CNVs (rare, may be artifacts)

**Rationale:**

Small CNVs (<30 kb) typically have lower Mendelian Precision due to:

* Lower signal-to-noise ratio
* Difficulty distinguishing from technical noise
* Higher breakpoint uncertainty
* Example: Analyse CNV from ≥30 kb

.. important::

   **Do not apply the same quality score threshold to all CNV sizes.** Small CNVs 
   require more stringent filtering than large CNVs. Always examine MP stratified 
   by size before setting filters.

Quality score filters
~~~~~~~~~~~~~~~~~~~~~

**Available metrics:**

* **Score** — Caller-specific quality score (higher = more confident)
* **SNP** — Number of SNP probes supporting the CNV (array data)
* **%overlap** — Reciprocal overlap percentage between CNVs from different algorithms
* **NbreAlgos** — Number of algorithms that detected the CNV (1, 2, 3, etc.)

**Rationale:**
Higher quality scores correlate with higher Mendelian Precision. The relationship is often non-linear, with MP plateauing at a certain threshold.

**Typical strategy:**

1. Plot MP versus quality score threshold
2. Observe the trade-off between quality and quantity
3. Select the quality score that meets your requirements

**Example:**

.. code-block:: text

   Deletions 50-100kb:
   
   Score ≥10  → MP = 60% (n = 800)
   Score ≥50  → MP = 75% (n = 500)
   Score ≥100 → MP = 85% (n = 300)
   Score ≥150 → MP = 91% (n = 150) 
   Score ≥200 → MP = 92% (n = 145)  
   
   Optimal threshold: 150 (trade-off between MP (quality) and CNV count (quantity))

**Caller concordance (if available):**

When merging CNV callsets from multiple algorithms, two metrics quantify caller agreement:

**1. NbreAlgos** — Number of algorithms detecting the CNV

* NbreAlgos = 1 → Detected by single algorithm only
* NbreAlgos = 2 → Detected by 2 algorithms
* NbreAlgos = 3 → Detected by 3 algorithms

**2. %overlap** — Reciprocal overlap percentage between detections

* %overlap = 0% → No overlap (NbreAlgos = 1)
* %overlap = 50% → 50% reciprocal overlap between algorithm calls
* %overlap = 100% → Perfect overlap between algorithm calls

**Relationship:**

* If **NbreAlgos = 1** → **%overlap = 0%** (no concordance possible)
* If **NbreAlgos ≥ 2** → **%overlap** can range from 0.1% to 100%

**Filtering strategy:**

Apply a minimum **%overlap** threshold to require caller concordance:

.. code-block:: text

   Filter                Meaning                              
   %overlap ≥ 0%        All CNVs (including single-caller)   
   %overlap ≥ 50%       ≥2 algos with ≥50% overlap           
   %overlap ≥ 70%       ≥2 algos with ≥70% overlap           
   %overlap ≥ 90%       ≥2 algos with ≥90% overlap           

**Note:** Filtering by **%overlap ≥ 50%** implicitly requires **NbreAlgos ≥ 2** 
(since single-caller CNVs have 0% overlap).

**Typical strategy:**

* **Balanced approach:** %overlap ≥ 50% (moderate concordance, retains more CNVs)
* **High precision:** %overlap ≥ 70% or ≥ 80% (strong concordance required)

**Rationale:**

CNVs detected independently by multiple algorithms with high reciprocal overlap 
are more likely to be genuine. Each caller has different sensitivities to 
artifacts, so concordance helps filter out caller-specific false positives.

.. tip::

   When working with merged callsets, **caller concordance (%overlap ≥ 50%)** 
   is often an effective filter. Apply this before optimizing other quality scores.

Problematic region filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Genomic regions prone to artifacts:**

* **Segmental duplications** — Highly similar sequences causing misalignment
* **Centromeres** — Repetitive, poorly mappable
* **Telomeres** — Repetitive, highly variable
* **HLA region** — Extreme polymorphism
* **Low mappability regions** — Reads cannot be uniquely placed

.. important::

   **Highly recommended:** Apply problematic region filters to all CNV datasets. 
   
   CNVs overlapping these regions have substantially lower Mendelian Precision 
   due to technical artifacts:
   
   * Read mismapping to paralogous sequences (segmental duplications)
   * Low coverage and poor mappability (centromeres, telomeres)
   * High genuine copy number variation (HLA region)

**Filter approach:**

* **Percent overlap threshold:** Exclude CNVs with >X% overlap with problematic regions
* **Binary filter:** Exclude any CNV overlapping problematic regions
* **Recommended strategy:** Apply 50% threshold (exclude CNVs with >50% overlap)

Transcript overlap filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Percent transcript overlap:**

* Exclude CNVs with low genic content (e.g., <10% overlap with transcripts)

**Use case:**

* When prioritizing functional variants

Gene-based filters
~~~~~~~~~~~~~~~~~~

**Exclusion lists:**

Upload a file of Ensembl Gene IDs to exclude **from MP calculation only**.

**Purpose:**

Based on published studies, you may identify genes known to be highly constrained and enriched for genuine *de novo* CNVs. Excluding these genes from MP calculation helps assess **technical precision** separately from biological *de novo* events.

**Examples:**

* Severe neurodevelopmental genes (*MECP2*, *SCN1A*, *CDKL5*)
* Haploinsufficient genes from disease databases (ClinGen, DDG2P)
* Genes with extreme constraint (pLI > 0.99)

**File format:**

Plain text file with one Ensembl Gene ID per line (no header):

.. code-block:: text

   ENSG00000169057
   ENSG00000198712
   ENSG00000130164

.. important::

   **Exclusion is for MP assessment only**
   
   These CNVs must be **retained** in your final dataset for downstream analyses 
   as they may represent pathogenic *de novo* variants.

**Workflow:**

1. Calculate MP (all CNVs) → e.g., 82%
2. Calculate MP (excluding gene list) → e.g., 90%
3. Difference (8%) estimates *de novo* contribution

**Gene constraint filters (LOEUF):**

Exclude CNVs affecting highly constrained genes (LOEUF < threshold, e.g., 0.6) **from MP calculation only**.

**Rationale:**

CNVs affecting constrained genes are enriched for genuine *de novo* events. These events reduce MP but are biologically valid. Excluding them from MP calculation allows you to:

1. **Assess technical precision** — MP without likely *de novo* events
2. **Estimate *de novo* rate** — Difference between filtered and unfiltered MP

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

**Typical workflow:**

1. Calculate MP (all CNVs) → e.g., 85% (technical + biological)
2. Calculate MP (LOEUF ≥ 0.6) → e.g., 92% (technical only)
3. Difference (92% - 85% = 7%) estimates genuine *de novo* contribution

MP Exploration filters
----------------------

**Purpose:**

Assess baseline MP and identify broad filtering needs.

**Available filters:**

* CNV size (slider)
* Minimum % transcript overlap (slider)
* Maximum % problematic regions overlap (slider)
* Gene exclusion list (upload)
* LOEUF threshold (slider)

**Workflow:**

1. Load preprocessed file
2. Choose transmission type (CNV-level or Gene-level)
3. Apply filters via sliders
4. Observe impact on MP (via plots and summary cards)
5. Proceed to Fine-tuning for detailed optimization

**Output:**

* MP by size range or MP by quality score stratified by size range
* Filtered table with Download option

Fine-tuning filters
-------------------

**Purpose:**

Systematically optimize MP.

**Available filters:**

* **CNV type:** DEL or DUP (analyzed separately)
* **Quality metrics:** Score, SNP, % overlap, nbre_algo
* **Operators:** ≥, ≤, = (combine multiple conditions)
* **Additional filters:** bp_overlap, LOEUF , t_Stop, t_Start

**Workflow:**

1. Apply filters
3. Compare "Before" vs "After" plots
4. Evaluate subset analyses (Genic only, Intergenic only, etc.)
5. Iterate until optimal threshold identified

**Output:**

* 4 comparative plots:
  * Before additional filters
  * After additional filters
  * After + Genic CNVs only
  * After + Intergenic CNVs only
* Downloadable tables for each plot

**Subset analyses:**

* **Genic CNVs only** — CNVs overlapping at least one gene
* **Intergenic CNVs only** — CNVs in non-genic regions
* **No excluded genes** — Exclude CNVs in user-provided gene list
* **No constrained genes (LOEUF < 1)** — Exclude CNVs in constrained genes

**Use cases:**

* Compare genic vs intergenic MP
* Assess impact of gene exclusion lists
* Evaluate technical MP (excluding constrained genes)

Optimization strategies
-----------------------

Plateau-based optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Approach:**

1. Plot MP versus quality score threshold
2. Identify plateau (where MP stops increasing)
3. Use lowest threshold at plateau

**Example:**

.. code-block:: text

   Threshold → MP      → CNVs retained
   ≥50       → 75%     → 1000
   ≥100      → 85%     → 600
   ≥150      → 92%     → 300  ← Plateau starts
   ≥200      → 92%     → 280  ← No further MP gain
   
   Optimal threshold: 150

**Rationale:**

Beyond the plateau, additional filtering removes genuine CNVs without improving MP.

Size-specific optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Approach:**

Optimize filters separately for each size range. Different size categories require different quality score thresholds.

.. important::

   **Critical:** Do not apply the same quality threshold to all CNV sizes. 
   This is a common mistake that leads to either:
   
   * Over-filtering large CNVs (losing genuine events)
   * Under-filtering small CNVs (retaining false positives)

**Example:**

.. code-block:: text

   Size range      Optimal Score threshold    MP after filtering
   1-30kb          ≥200                       85%
   30-50kb         ≥150                       90%
   50-100kb        ≥100                       92%
   100-200kb       ≥50                        94%
   >200kb          ≥20                        95%

**Rationale:**

Small CNVs have lower signal-to-noise ratio and require more stringent filtering to achieve comparable MP. Large CNVs are inherently more reliable and can pass with lower quality scores.

**Workflow:**

1. Stratify MP by size range
2. For each size range, plot MP vs quality threshold
3. Identify plateau for each size range
4. Apply size-specific thresholds

Type-specific optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Approach:**

Optimize filters separately for deletions (DEL) and duplications (DUP).

.. tip::

   Always analyze DEL and DUP separately. They have different baseline MP values 
   and respond differently to filtering.

**Typical patterns:**

* **Deletions:** Higher baseline MP (signal easier to detect)
* **Duplications:** Lower baseline MP (signal harder to detect, more ambiguous)

**Implications for filtering:**

* Deletions may achieve high MP (≥90%) with moderate filtering
* Duplications may require more aggressive filtering to reach similar MP

**Example:**

.. code-block:: text

   CNV type    Size range    Optimal threshold    MP after filtering
   DEL         50-100kb      Score ≥100           92%
   DUP         50-100kb      Score ≥150           88%

**Rationale:**

Duplications are technically harder to detect than deletions. CNV callers typically have:

* Higher sensitivity for deletions (easier to detect copy number loss)
* Lower sensitivity for duplications (copy number gain harder to distinguish from noise)

This difference in detection difficulty translates to different optimal filtering thresholds.

Gene constraint consideration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Approach:**

When evaluating filtering strategies, separate technical precision from biological *de novo* rate using LOEUF.

**Workflow:**

1. Apply your candidate filter strategy
2. Calculate MP (all CNVs)
3. Calculate MP (excluding LOEUF < 0.6)
4. If the difference is large (>10%), many non-inherited CNVs may be genuine *de novo* events

**Example:**

.. code-block:: text

   Filter: Score ≥100, Size ≥30kb
   
   MP (all CNVs) = 85%
   MP (LOEUF ≥ 0.6) = 92%
   
   Interpretation:
   - Technical precision: 92% (good)
   - Estimated *de novo* rate: 7% (reasonable)
   - Filter strategy is appropriate

vs.

.. code-block:: text

   Filter: Score ≥50, Size ≥10kb
   
   MP (all CNVs) = 70%
   MP (LOEUF ≥ 0.6) = 72%
   
   Interpretation:
   - Technical precision: 72% (poor)
   - Estimated *de novo* rate: 2%
   - High false positive rate → More aggressive filtering needed

**Rationale:**

This approach helps you distinguish:

* Low MP due to technical false positives (requires filtering)
* Low MP due to genuine *de novo* events (biologically expected)

Balancing precision and sensitivity
------------------------------------

**Trade-off:**

* **More aggressive filtering** → Higher MP, fewer CNVs
* **Lenient filtering** → Lower MP, more CNVs

**Considerations:**

* **For discovery studies:** Prioritize sensitivity (lenient filtering, MP ≥70%)
* **For clinical validation:** Prioritize precision (aggressive filtering, MP ≥90%)
* **For method comparison:** Use consistent filters across methods

**Recommended approach:**

1. Start with minimal filters to assess baseline quality
2. Identify optimal thresholds via Fine-tuning 
3. Apply filters that achieve MP ≥85% while retaining sufficient CNVs
4. For clinical applications, target MP ≥90%

.. admonition:: Best practices

  1. **Always stratify by size** before filtering — Do not apply uniform thresholds
  2. **Always calculate MP separately for DEL and DUP** — They have different quality profiles
  3. **Apply problematic region filters** — Highly recommended for all datasets
  4. **Consider gene constraint (LOEUF)** — Distinguish technical FP from biological *de novo*
  5.  **Caller concordance (if available):** Require ≥2 algorithms with ≥50% reciprocal overlap 
  6. **Balance quality and quantity** — Visualize the trade-off between MP (quality) and CNV count (quantity) to make informed filtering decisions based on your study goals

See also
--------

* :doc:`mendelian_precision` — Understanding the MP metric
* :doc:`preprocessing` — Pre-filtering annotation and inheritance calculation
* :doc:`outputs` — Visualizing filtering impact and downloading filtered tables
