Inheritance status
==================

MCNV2 classifies each offspring CNV as **inherited (True)** or **non-inherited (False)** based on parental CNV data.

This page explains the two inheritance modes available in MCNV2.

Overview
--------

Inheritance status is calculated using two approaches:

1. **CNV-level matching** (coordinate-based) → Column: **Transmitted_CNV**
2. **Gene-level matching** (gene-based) → Column: **Transmitted_gene**

Both approaches are computed during the Preprocessing step and stored in the inheritance status table.

CNV-level matching (Transmitted_CNV)
-------------------------------------

**Definition:**

A child CNV is considered **inherited (True)** if it has reciprocal overlap ≥ threshold with at least one parental CNV (same type and chromosome).

**Algorithm:**

For each child CNV:

1. Find all parental CNVs (father and mother) on the same chromosome with the same type (DEL or DUP)
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

* Sensitive to breakpoints
* May miss true inheritances if parental CNV has slightly different boundaries

Gene-level matching (Transmitted_gene)
---------------------------------------

**Definition:**

A child CNV is considered **inherited (True)** if at least one affected gene is also affected in at least one parent.

**Algorithm:**

For each child CNV-gene pair:

1. Identify all genes overlapping the child CNV
2. For each gene, check if it is also affected in the father or mother
3. If yes → **True** (inherited)
4. If no → **False** (non-inherited)
5. If CNV is intergenic (no gene overlap) → **"intergenic"**

**Example:**

.. code-block:: text

   Child CNV:   chr1:1,000,000-1,100,000 → Overlaps genes A, B, C
   Father CNV:  chr1:1,050,000-1,150,000 → Overlaps genes B, C, D
   Mother:      No CNVs affecting genes A, B, C

   Gene A: Not in father → Transmitted_gene = False
   Gene B: In father → Transmitted_gene = True
   Gene C: In father → Transmitted_gene = True

**Advantages:**

* Robust to breakpoint variability 
* Focuses on functional impact (affected genes)
* Useful for evaluating gene-based filters

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
| **Threshold**             | User-defined (e.g., 50%)| Not applicable         |
+---------------------------+-------------------------+------------------------+
| **Applicable to**         | All CNVs                | Genic CNVs only        |
+---------------------------+-------------------------+------------------------+
| **Robustness to**         | Low (sensitive to       | High (ignores          |
| **breakpoint shifts**     | boundaries)             | boundaries)            |
+---------------------------+-------------------------+------------------------+


See also
--------

* :doc:`preprocessing` — How inheritance status is calculated during preprocessing
* :doc:`mendelian_precision` — How inheritance status is used to compute MP
* :doc:`filtering` — How to filter CNVs based on inheritance status
