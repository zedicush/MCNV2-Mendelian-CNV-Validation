#!/usr/bin/env python3
import argparse
import polars as pl


"""
Annotate CNVs with Transmitted Gene Information

This script reads a filtered pedigree file and annotated CNVs, identifies parent and child CNVs,
and flags which disrupted genes are transmitted from parents to children.

Usage example:
    python get_transmitted_genes.py \
        --pedigree data/filtered_pedigree.tsv \
        --anno_cnv_file data/annotated_cnvs.tsv \
        --output cnv_anno_with_transmitted_gene.tsv

Dependencies:
    - Python 3.8+
    - Polars 1.27+
    
Author:
    Benjamin Clark
"""
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Read pedigree and CNV annotation files")
    parser.add_argument(
        "--pedigree",
        type=str,
        required=True,
        help="Path to the filtered pedigree TSV file"
    )
    parser.add_argument(
        "--anno_cnv_file",
        type=str,
        required=True,
        help="Path to the annotated CNVs TSV file"
    )

    parser.add_argument(
        "--output",
        type=str,
        required=False,
        default="cnv_anno_with_transmitted_gene.tsv",
        help="Output filename"
    )

    args = parser.parse_args()

    # Read CSVs
    ped = pl.scan_csv(args.pedigree, separator='\t')
    anno_cnv = pl.scan_csv(args.anno_cnv_file, separator='\t').drop_nulls()


    #Defining parent CNVs
    parents_cnv = pl.concat([anno_cnv.join(ped, left_on='SampleID', right_on=["FatherID"], how='semi'),
                             anno_cnv.join(ped, left_on='SampleID', right_on=["MotherID"], how='semi')])
    #Defining child CNVs
    child_cnv = anno_cnv.join(ped, how='semi', on='SampleID')

    #Defining transmitted disrupted genes
    inherited_cnv = pl.concat([(child_cnv.join(parents_cnv, how='semi', on=['SampleID', "gene_name", "TYPE"])  #inherited are found in the parent table
                              .with_columns(pl.lit(True).alias("transmittedGene"))),
                        child_cnv.join(parents_cnv, how='anti', on=['SampleID', "gene_name", "TYPE"])          #not inherited are not found in the parent table
                              .with_columns(pl.lit(False).alias("transmittedGene"))])
    #writing out
    inherited_cnv.sink_csv(f"{args.output}", separator='\t')

if __name__ == "__main__":
    main()
