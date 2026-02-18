Outputs
=======

MCNV2 generates two types of outputs: **visualizations** (plots) for exploring Mendelian Precision patterns, and **downloadable tables** for downstream analyses.

Visualizations
--------------

MCNV2 provides interactive plots to assess CNV quality and identify optimal filtering strategies.

MP by size range (bar plots)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Description:**

Bar plots showing Mendelian Precision for each CNV size range.

**Axes:**

* **X-axis:** Size ranges (1-30kb, 30-50kb, 50-100kb, 100-200kb, 200-500kb, 500kb-1Mb, >1Mb)
* **Y-axis:** Mendelian Precision (%)

**Data displayed:**

* Separate plots for deletions (DEL) and duplications (DUP)
* Bar height represents MP
* Number on top of each bar shows CNV count (n)

**Use case:**

Identify size ranges with low MP that may benefit from filtering or exclusion.

MP versus quality score (line plots)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Description:**

Line plots showing Mendelian Precision as a function of quality score threshold stratified by size range.

**Axes:**

* **X-axis:** Quality score threshold (e.g., Score ≥ 0, 10, 20, ..., 200)
* **Y-axis:** Mendelian Precision (%)

**Data displayed:**

* One line per size range (color-coded as above)
* Each point represents MP when applying quality score threshold ≥ X
* Separate plots for deletions (DEL) and duplications (DUP)

**Legend:**

Right side of plot shows:

* **n** — Total CNV count for each size range
* **CNV Length** — Size range label with marker symbol
* Marker symbols: circles, triangles, plus signs, crosses, diamonds, etc.

**Interpretation:**

* **Plateau pattern:** MP increases rapidly, then plateaus
* **Optimal threshold:** Lowest value where plateau is reached
* **Divergent curves:** Different size ranges require different thresholds

**Interactive features:**

* **Hover:** Tooltip shows threshold value, MP, size range, and CNV count
* **Click legend:** Show/hide specific size range curves
* **Zoom:** Select area to zoom in
* **Download:** Export as PNG, PDF, or SVG

**Use case:**

Identify optimal quality score thresholds for each size range.

Plot interactivity
~~~~~~~~~~~~~~~~~~

**Hover tooltips:**

When hovering over data points, a tooltip displays:

* **threshold:** Quality score value
* **MP:** Mendelian Precision at this threshold
* **Size_Range:** CNV size category
* **n:** Number of CNVs remaining after applying filter

**Multiple simultaneous tooltips:**

Hovering near multiple curves shows tooltips for all size ranges at that threshold, enabling direct comparison.

**Legend interaction:**

Click on size range labels in the legend to show/hide curves. Useful for focusing on specific size categories.

**Zoom and pan:**

* **Box select:** Click and drag to zoom into a region
* **Reset:** Double-click to reset zoom
* **Pan:** Shift + drag to move view

**Download options:**

Click camera icon to download plot as:

* PNG (raster image)

Comparative plots (Fine-tuning)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Before vs After:**

Fine-tuning displays up to 4 plots simultaneously:

1. **Before additional filters** — Baseline MP (from MP Exploration)
2. **After additional filters** — MP after applying Fine-tuning filters
3. **After + subset 1** — MP after filters + subset (e.g., Genic only)
4. **After + subset 2** — MP after filters + subset (e.g., Intergenic only)

**Purpose:**

* Assess impact of additional filters
* Compare genic vs intergenic MP
* Evaluate subset-specific patterns

**Zoom modal:**

Click **+** button on any plot to view it enlarged in a modal window with:

* Full-screen visualization
* Enhanced interactivity
* Download table button (green)
* Close button

Downloadable tables
-------------------

MCNV2 provides downloadable tables at multiple stages of analysis.

MP Exploration table
~~~~~~~~~~~~~~~~~~~~

**Location:**

MP Exploration interface → "Filtered table" tab → Download CSV button

**Content:**

All CNVs passing the applied filters (size, transcript overlap, problematic regions, LOEUF).

**Columns:**

* **Original CNV columns:** CHR, START, STOP, TYPE, SAMPLE_ID, quality scores
* **Annotation columns:** GeneName, GeneID, Transcript, LOEUF, problematic_region_overlap
* **Inheritance columns:** Transmitted_CNV (True/False), Transmitted_gene (True/False/intergenic)
* **Additional columns:** TrioKey, family_statue, cnv_id, Size_Range, transmission

**Row structure:**

* One row per CNV-gene pair (CNVs overlapping multiple genes have multiple rows)
* Intergenic CNVs have empty gene columns

**Use case:**

* Downstream statistical analyses
* Manual inspection of specific CNVs
* Integration with other genomic datasets

Export workflow
---------------

MP Exploration
~~~~~~~~~~~~~~

**Step 1:** Apply filters using sliders and dropdowns

**Step 2:** Click "Apply filters" button

**Step 3:** Navigate to "Filtered table" tab

**Step 4:** Click "Download CSV" button

**Step 5:** Save file to desired location

**Filename format:** `cnvs_filtered_<timestamp>.csv`

Fine-tuning
~~~~~~~~~~~

**Step 1:** Apply advanced filters using dropdowns and value fields

**Step 2:** Click "Apply filters" button

**Step 3:** Click **+** button on desired plot (Before, After, or subset)

**Step 4:** In enlarged modal, click green "Download table" button

**Step 5:** Save file to desired location

**Filename format:** `cnvs_finetuned_<condition>_<timestamp>.csv`

See also
--------

* :doc:`mendelian_precision` — Understanding MP metrics
* :doc:`filtering` — Filtering strategies
* :doc:`preprocessing` — Pre-export annotation
