Plots and visualizations
========================

MCNV2 generates two types of Mendelian Precision (MP) visualizations to help identify high-quality CNV calls and optimal filtering thresholds.

MP by size range (bar plot)
----------------------------

**Type:** Bar plot  
**X-axis:** CNV size categories  
**Y-axis:** Mendelian Precision (%)

**Size categories:**

* <30 kb
* 30–50 kb
* 50–100 kb
* 100–500 kb
* 500 kb–1 Mb
* >1 Mb

**Purpose:**

This plot shows how Mendelian Precision varies with CNV size. Typically:

* Small CNVs (<100 kb) have lower baseline MP (~35-72% for deletions)
* Large CNVs (>1 Mb) have high MP (>95%)

**Example interpretation:**

.. code-block:: text

   Size range       MP (deletions)
   <30 kb           35%
   30-50 kb         45%
   50-100 kb        60%
   100-500 kb       82%
   500kb-1Mb        91%
   >1 Mb            97%

**Key insight:** Larger CNVs are generally more reliable (higher MP) because they are easier to detect and less affected by technical noise.

.. note::

   This visualization is **always available**, regardless of whether quality scores are provided in the CNV file.

MP versus quality score (line plot)
------------------------------------

**Type:** Line plot (one curve per size range)  
**X-axis:** Quality score (caller-specific)  
**Y-axis:** Mendelian Precision (%)

**Purpose:**

This plot identifies the **optimal quality score threshold** where MP plateaus, indicating high-confidence calls.

**Requirements:**

.. important::

   This visualization requires **quality scores** (column 6) in the CNV file.
   
   Without quality scores, this plot cannot be generated.

**Typical pattern:**

.. code-block:: text

   Quality Score    MP (deletions 30-100kb)
   0                35%
   10               42%
   20               58%
   30               72%  ← MP plateaus
   40               73%
   50               73%

**Key insight:** For this example, setting a quality score threshold of ≥30 captures most high-confidence calls without sacrificing too much sensitivity.

.. tip::

   Use this plot to determine caller-specific quality thresholds. Different CNV callers (PennCNV, GATK-gCNV, CNVnator) have different quality score scales.

Combining both visualizations
------------------------------

**Recommended workflow:**

1. **Start with size-stratified bar plots** to understand baseline MP across size ranges
2. **Use MP vs quality score curves** to identify optimal thresholds for each size category
3. **Apply combined filters** (quality score + size + other filters) to maximize MP while maintaining yield

**Example decision process:**

.. code-block:: text

   Small deletions (30-100 kb):
   - Baseline MP: 60%
   - MP at QS≥30: 72%
   - Decision: Apply QS≥30 filter → +12% MP gain

   Large deletions (>1 Mb):
   - Baseline MP: 97%
   - MP at QS≥30: 97%
   - Decision: No need for strict QS filter (already high quality)

Interactive exploration (Shiny app)
------------------------------------

The Shiny application allows **real-time filtering** and automatic plot updates:

* Adjust quality score thresholds with sliders
* Toggle filters (caller concordance, LOEUF, problematic regions)
* See MP recalculated instantly
* Compare filtering strategies side-by-side

**Workflow:**

1. Upload CNV and pedigree files
2. View baseline MP (bar plots + curves)
3. Apply filters interactively
4. Observe MP changes in real-time
5. Export final filtered dataset and publication-ready figures

.. tip::

   The Shiny interface is ideal for exploratory analysis and threshold optimization. Once optimal parameters are identified, use the CLI for batch processing of additional cohorts.

Exporting figures
-----------------

All plots can be exported in publication-ready formats:

* **PNG** (high resolution, 300 dpi)
* **PDF** (vector graphics, scalable)
* **SVG** (vector graphics, editable)

**From Shiny app:** Click "Download plot" button  
**From CLI:** Specify output format in command arguments

See :doc:`exports` for detailed export options.
