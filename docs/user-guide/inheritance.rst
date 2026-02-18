Inheritance status
==================

MCNV2 classifies each offspring CNV as **inherited (True)** or **non-inherited (False)** based on parental CNV data.

This page provides a comprehensive explanation of the two inheritance modes available in MCNV2.

Overview
--------

Inheritance status is calculated using two approaches:

1. **CNV-level matching** (coordinate-based) → Column: **Transmitted_CNV**
2. **Gene-level matching** (gene-based) → Column: **Transmitted_gene**

Both approaches are computed during the Preprocessing step and stored in the inheritance status table.

CNV-level matching (Transmitted_CNV)
-------------------------------------
.. note:: 
   **Definition:**

   A child CNV is considered **inherited (True)** if it has reciprocal overlap ≥ threshold with at least one parental CNV of the **same type (DEL or DUP)** on the same chromosome.

**Algorithm:**

For each child CNV:

1. Find all parental CNVs (father and mother) on the same chromosome with the **same type (DEL or DUP)**
2. Calculate reciprocal overlap with each parental CNV:

   .. math::

      \text{Reciprocal overlap} = \frac{\text{Overlap length}}{\min(\text{Child CNV length}, \text{Parent CNV length})}

3. If reciprocal overlap ≥ threshold (e.g., 50%) with at least one parent → **True** (inherited)
4. Otherwise → **False** (non-inherited, candidate *de novo*)

**Example:**

.. code-block:: text

   Child CNV:   chr1:1,000,000-1,100,000 (100 kb deletion)
   Father CNV:  chr1:1,020,000-1,120,000 (100 kb deletion)
   Mother:      No overlapping CNV

   Overlap length: 80 kb (1,020,000 to 1,100,000)
   Min length: 100 kb
   Reciprocal overlap: 80 / 100 = 0.80 (80%)

   Threshold: 0.5 (50%)
   Result: 80% ≥ 50% → Transmitted_CNV = True (inherited from father)

**Advantages:**

* Simple and intuitive
* Works for all CNVs (genic and intergenic)
* Standard approach in CNV literature

**Limitations:**

* Requires defining an overlap threshold (no consensus on optimal value)
* Sensitive to breakpoint differences — may miss true inheritances when child and parent CNVs have slightly different boundaries
* Sensitive to CNV fragmentation — a CNV detected as multiple fragments may not be recognized as inherited even if it represents the same parental CNV (e.g., child fragments: 150kb + 200kb + 50kb vs parent: single 400kb call)

Gene-level matching (Transmitted_gene)
---------------------------------------
.. note:: 
   **Definition:**

   A child CNV is considered **inherited (True)** if at least one affected gene is also affected by the **same CNV type (DEL or DUP)** in at least one parent.

**Algorithm:**

For each child CNV-gene pair:

1. Identify all genes overlapping the child CNV
2. For each gene, check if it is also affected by the **same CNV type** in the father or mother
3. If yes → **True** (inherited)
4. If no → **False** (non-inherited)
5. If CNV is intergenic (no gene overlap) → **"intergenic"**

**Example 1: Same type (inherited):**

.. code-block:: text

   Child CNV:   chr1:1,000,000-1,100,000 (DEL, overlaps genes A, B, C)
   Father CNV:  chr1:1,050,000-1,150,000 (DEL, overlaps genes B, C, D)
   Mother:      No CNVs affecting genes A, B, C

   Gene A: Not in father (DEL) → Transmitted_gene = False
   Gene B: In father (DEL) → Transmitted_gene = True
   Gene C: In father (DEL) → Transmitted_gene = True

**Example 2: Different type (not inherited):**

.. code-block:: text

   Child CNV:   chr1:1,000,000-1,100,000 (DEL, overlaps gene BRCA1)
   Father CNV:  chr1:1,020,000-1,150,000 (DUP, overlaps gene BRCA1)
   Mother:      No CNVs

   Gene BRCA1 is deleted in child but duplicated in father
   → Different CNV types
   → Transmitted_gene = False

**Advantages:**

* Robust to breakpoint variability and CNV fragmentation
* Focuses on functional impact (affected genes)
* Useful for evaluating gene-based filters
* No threshold parameter required

**Limitations:**

* Only applicable to genic CNVs
* Intergenic CNVs are flagged as "intergenic" (not True/False)

Comparison of the two approaches
---------------------------------

+---------------------------+-------------------------+------------------------+
| Aspect                    | Transmitted_CNV         | Transmitted_gene       |
+===========================+=========================+========================+
| **Method**                | Coordinate overlap      | Shared genes           |
+---------------------------+-------------------------+------------------------+
| **Threshold**             | User-defined (e.g., 50%)| Not required           |
+---------------------------+-------------------------+------------------------+
| **Applicable to**         | All CNVs                | Genic CNVs only        |
+---------------------------+-------------------------+------------------------+
| **Robustness to**         | Low (sensitive to       | High (ignores          |
| **breakpoint shifts**     | boundaries)             | boundaries)            |
+---------------------------+-------------------------+------------------------+
| **Robustness to**         | Low (fragments may      | High (focuses on       |
| **CNV fragmentation**     | fail threshold)         | shared genes)          |
+---------------------------+-------------------------+------------------------+

.. tip::
  **Recommendation:**

  Use both approaches to get a comprehensive view of CNV quality. The difference between CNV-level and gene-level MP highlights the impact of breakpoint variability and CNV fragmentation in your dataset.

Example use cases
-----------------

**Case 1: Inherited by both methods**

.. code-block:: text

   Child:  chr1:1,000,000-1,100,000 (DEL, overlaps gene X)
   Father: chr1:1,000,000-1,100,000 (DEL, overlaps gene X)

   Transmitted_CNV: True (100% overlap)
   Transmitted_gene: True (gene X deleted in father)

**Interpretation:** High-confidence inherited CNV.

**Case 2: Inherited by gene, not by CNV (breakpoint shift)**

.. code-block:: text

   Child:  chr1:1,000,000-1,100,000 (DEL, overlaps gene X)
   Father: chr1:1,200,000-1,300,000 (DEL, overlaps gene X)

   Transmitted_CNV: False (no coordinate overlap)
   Transmitted_gene: True (gene X deleted in father)

**Interpretation:** Different breakpoints, but same gene affected. Likely inherited despite coordinate mismatch.

**Case 3: Inherited by gene, not by CNV (fragmentation)**

.. code-block:: text

   Child:  chr1:1,000,000-1,150,000 (DEL, 150kb, overlaps gene X)
           chr1:1,150,000-1,350,000 (DEL, 200kb, overlaps gene X)
           chr1:1,350,000-1,400,000 (DEL, 50kb, no gene)
   
   Father: chr1:1,000,000-1,300,000 (DEL, 300kb, overlaps gene X)

   Transmitted_CNV: Fragment 1 (100% overlap) → True
                    Fragment 2 (75% overlap) → True if threshold ≤75%, False if >75%
                    Fragment 3 (0% overlap) → False
   
   Transmitted_gene: All fragments overlapping gene X → True

**Interpretation:** CNV fragmentation causes some fragments to fail coordinate matching, but gene-level matching correctly identifies inheritance.

**Case 4: Non-inherited (candidate *de novo*)**

.. code-block:: text

   Child:  chr1:1,000,000-1,100,000 (DEL, overlaps gene X)
   Father: No CNV affecting gene X
   Mother: No CNV affecting gene X

   Transmitted_CNV: False
   Transmitted_gene: False

**Interpretation:** Candidate *de novo* CNV or technical artefact.

**Case 5: Intergenic CNV**

.. code-block:: text

   Child:  chr1:10,000,000-10,050,000 (DEL, intergenic)
   Father: chr1:10,000,000-10,050,000 (DEL, intergenic)

   Transmitted_CNV: True (100% overlap)
   Transmitted_gene: "intergenic"

**Interpretation:** Inherited based on coordinates. Gene-based method not applicable.

**Case 6: Different CNV type**

.. code-block:: text

   Child:  chr1:1,000,000-1,100,000 (DEL, overlaps gene X)
   Father: chr1:1,000,000-1,100,000 (DUP, overlaps gene X)

   Transmitted_CNV: False (different types)
   Transmitted_gene: False (different types)

**Interpretation:** Same coordinates but opposite CNV types. Not inherited.


See also
--------

* :doc:`preprocessing` — How inheritance status is calculated during preprocessing
* :doc:`mendelian_precision` — How inheritance status is used to compute MP
* :doc:`filtering` — How to filter CNVs based on inheritance status
