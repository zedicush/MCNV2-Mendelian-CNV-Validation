Definitions
===========

| **Copy Number Variation CNV**: CNV are deletions or duplications of genomic regions greater than or equal to 1000 nucleotides and can encompass one or more genes. 
| **Trio**: offspring + mother + father.
| **Inherited**: present in offspring and ≥1 parent.
| **Mendelian error**: present in offspring, absent in both parents.
| **Mendelian Precision (MP)**: proportion of inherited CNVs among all offspring CNVs.
| **Reciprocal overlap**: overlap(child,parent)/len(child) and overlap(child,parent)/len(parent).
| **LOEUF**: gnomAD constraint score; lower values indicate stronger constraint.

--------------

MP formula
----------
.. math::
   :label: mp-formula

   MP = \frac{I}{N} 

Where:
- I = Inherited CNVs
- N = Total CNVs  

