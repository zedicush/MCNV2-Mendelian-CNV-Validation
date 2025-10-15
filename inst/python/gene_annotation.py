#!/usr/bin/env python3
import argparse
import polars as pl
import subprocess
import tempfile
import shutil
import sys

"""
CNV–Gene Intersection (Polars + bedtools) avec ordre de sortie figé.

Entrées
- --cnv : TSV. Les 5 premières colonnes (quel que soit leur nom) sont mappées positionnellement à :
          Chr, Start, Stop, Type, SampleID. Les colonnes suivantes sont conservées dans leur ordre.
- --gene_resource : ressource gène/GTFlike utilisée en -b de bedtools intersect (le script suppose
                    la même structure que celle avec laquelle vous avez déjà testé).
- --prob_regions : TSV/BED avec/ou sans header ; si GenomeVersion existe, filtrage par alias.
- --genome_version : alias parmi GRCh38, GRCh37, hg19, hg38, 19, 38.

Sortie (ordre de colonnes figé)
1) Chr, Start, Stop, Type, SampleID
2) Colonnes CNV supplémentaires (ordre d’origine)
3) Colonnes gènes/pRegions “fixes” (dans cet ordre) :
   gene_type, transcript, gene_name, bp_overlap, gene_id, LOEUF, t_Start, t_Stop
4) Autres colonnes issues du fichier gène (ordre d’apparition), en excluant t_End
5) cnv_problematic_region_overlap
"""

# ----- Helpers -----

def _ci_lookup(cols, wanted):
    w = wanted.lower()
    for c in cols:
        if c.lower() == w:
            return c
    return None

def _normalize_gv_arg(gv: str) -> str:
    s = (gv or "").strip().lower()
    if "grch38" in s or "hg38" in s or s == "38":
        return "38"
    if "grch37" in s or "hg19" in s or s in {"19", "37"}:
        return "37"
    digits = "".join(ch for ch in s if ch.isdigit())
    if digits == "38":
        return "38"
    if digits in {"37", "19"}:
        return "37"
    return "38"

def _pattern_for_gv(norm: str) -> str:
    return r"(grch38|hg38|38)" if norm == "38" else r"(grch37|hg19|37|19)"

def _load_cnv_positional(path: str) -> pl.LazyFrame:
    """
    Lit le CNV TSV. Quelle que soit l’entête :
    - mappe positionnellement les 5 premières colonnes -> Chr, Start, Stop, Type, SampleID
    - conserve toutes les colonnes suivantes dans leur ordre.
    """
    try:
        lf = pl.scan_csv(path, separator="\t", has_header=True)
        cols = lf.collect_schema().names()
        if len(cols) < 5:
            raise ValueError
        c1, c2, c3, c4, c5 = cols[0], cols[1], cols[2], cols[3], cols[4]
    except Exception:
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        cols = lf.collect_schema().names()
        if len(cols) < 5:
            raise ValueError("CNV file must have at least 5 columns (Chr, Start, Stop, Type, SampleID).")
        c1, c2, c3, c4, c5 = cols[0], cols[1], cols[2], cols[3], cols[4]

    cnv_extras = cols[5:]

    lf = (
        lf.select(
            [
                pl.col(c1).alias("Chr"),
                pl.col(c2).cast(pl.Int64, strict=False).alias("Start"),
                pl.col(c3).cast(pl.Int64, strict=False).alias("Stop"),
                pl.col(c4).cast(pl.Utf8,  strict=False).alias("Type"),
                pl.col(c5).cast(pl.Utf8,  strict=False).alias("SampleID"),
            ] + [pl.col(c) for c in cnv_extras]
        )
        .drop_nulls(["Start", "Stop"])
    )
    return lf, cnv_extras  # on renvoie aussi la liste des extras pour l’ordre final

def _load_prob_regions_flex(path: str, genome_version: str) -> pl.LazyFrame:
    """
    Lit prob_regions (TSV/BED) ; si header et colonne GenomeVersion → filtre par alias.
    Normalise en Chr, Start, Stop.
    """
    norm = _normalize_gv_arg(genome_version)
    pattern = _pattern_for_gv(norm)

    try:
        lf = pl.scan_csv(path, separator="\t", has_header=True)
        cols = lf.collect_schema().names()
        if len(cols) < 3:
            raise ValueError

        chr_col   = _ci_lookup(cols, "Chr")   or cols[0]
        start_col = _ci_lookup(cols, "Start") or cols[1]
        stop_col  = _ci_lookup(cols, "Stop") or _ci_lookup(cols, "End") or cols[2]

        gv_col = _ci_lookup(cols, "GenomeVersion")
        if gv_col:
            lf = lf.filter(pl.col(gv_col).cast(pl.Utf8).str.to_lowercase().str.contains(pattern))
        else:
            print("[WARN] 'GenomeVersion' not found in prob_regions. No filtering applied.", file=sys.stderr)

        lf = (
            lf.select([
                pl.col(chr_col).alias("Chr"),
                pl.col(start_col).cast(pl.Int64, strict=False).alias("Start"),
                pl.col(stop_col).cast(pl.Int64, strict=False).alias("Stop"),
            ])
            .drop_nulls(["Start", "Stop"])
        )
        return lf
    except Exception:
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        return (
            lf.select([
                pl.col("column_1").alias("Chr"),
                pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
            ])
            .drop_nulls(["Start", "Stop"])
        )

def _finalize_order(
    joined_df: pl.DataFrame,
    cnv_extras: list[str],
    gene_schema_cols: list[str]
) -> pl.DataFrame:
    """
    Fixe l’ordre final des colonnes :
      Front CNV (5) -> extras CNV -> gene fixed -> autres colonnes gènes (sans t_End) -> cnv_problematic_region_overlap
    """
    front = ["Chr", "Start", "Stop", "Type", "SampleID"]

    # colonnes gène/pRegions fixées, dans cet ordre exact
    gene_fixed = [
        "gene_type", "transcript", "gene_name", "bp_overlap", "gene_id", "LOEUF", "t_Start", "t_Stop"
    ]

    exist = set(joined_df.columns)

    # colonnes CNV extra (uniquement celles effectivement présentes après les joins)
    cnv_extra_kept = [c for c in cnv_extras if c in exist]

    # Toutes les colonnes provenant (potentiellement) du fichier gènes qu’on a vues dans le schéma
    # On supprime ce qui est déjà dans gene_fixed + front + extras, et on exclut explicitement t_End
    already = set(front) | set(cnv_extra_kept) | set(gene_fixed)
    gene_all_present = [c for c in gene_schema_cols if c in exist]
    gene_other = [c for c in gene_all_present if c not in already and c != "t_End"]

    # ‘cnv_problematic_region_overlap’ en dernier (si présent)
    tail = ["cnv_problematic_region_overlap"] if "cnv_problematic_region_overlap" in exist else []

    final_cols = [c for c in front if c in exist] \
               + cnv_extra_kept \
               + [c for c in gene_fixed if c in exist] \
               + gene_other \
               + tail

    return joined_df.select(final_cols)

# ----- Main pipeline -----

def main():
    ap = argparse.ArgumentParser(
        description="Build CNV–gene intersection DB with fixed output column order."
    )
    ap.add_argument("--cnv", required=True, help="CNV TSV. First 5 cols mapped to Chr, Start, Stop, Type, SampleID.")
    ap.add_argument("--gene_resource", required=True, help="Gene resource TSV/BED used with bedtools -b.")
    ap.add_argument("--prob_regions", required=True, help="Problematic regions TSV/BED.")
    ap.add_argument("--genome_version", default="GRCh38", help="Genome alias: GRCh38/GRCh37/hg19/hg38/19/38.")
    ap.add_argument("--bedtools_path", default="bedtools", help="Absolute path to bedtools executable")
    ap.add_argument("--out", default="data/cnv_geneDB.tsv", help="Output TSV path.")
    args = ap.parse_args()

    # CNV (positionnel) + prob_regions
    cnvLF, cnv_extras = _load_cnv_positional(args.cnv)
    pRegions = _load_prob_regions_flex(args.prob_regions, args.genome_version)

    tmpdir = tempfile.mkdtemp()
    try:
        # Écrit des BED 3-col pour bedtools
        cnv_bed = f"{tmpdir}/tmp_cnvs.bed"
        cnvLF.select(["Chr", "Start", "Stop"]).unique().sink_csv(cnv_bed, separator="\t", include_header=False)

        preg_bed = f"{tmpdir}/tmp_pRegions.bed"
        pRegions.select(["Chr", "Start", "Stop"]).sink_csv(preg_bed, separator="\t", include_header=False)

        # bedtools: CNV x Gene
        inter_bed = f"{tmpdir}/tmp_intersect.bed"
        subprocess.run(
            f"bedtools intersect -a {cnv_bed} -b {args.gene_resource} -F 0.1 -wao > {inter_bed}",
            shell=True, check=True
        )

        # bedtools: CNV x prob_regions
        preg_inter_bed = f"{tmpdir}/tmp_pRegion_intersect.bed"
        subprocess.run(
            f"bedtools intersect -a {cnv_bed} -b {preg_bed} -wao > {preg_inter_bed}",
            shell=True, check=True
        )

        # Parse intersection CNV x Gene (schéma correspondant à ta ressource actuelle)
        # Colonnes CNV: col1..3 ; Colonnes GENE: ici on suppose attrs=col7, puis colonnes “typiques” comme déjà validé
        gene_lf = (
            pl.scan_csv(inter_bed, separator="\t", has_header=False)
              .select([
                  pl.col("column_1").alias("Chr"),
                  pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                  pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
                  pl.col("column_5").cast(pl.Int64, strict=False).alias("t_Start"),
                  pl.col("column_6").cast(pl.Int64, strict=False).alias("t_End"),
                  pl.col("column_7").cast(pl.Utf8,   strict=False).alias("attrs"),
                  pl.col("column_8").cast(pl.Utf8,   strict=False).alias("gene_type"),
                  pl.col("column_9").cast(pl.Utf8,   strict=False).alias("transcript"),
                  pl.col("column_10").cast(pl.Utf8,  strict=False).alias("gene_name"),
                  pl.col("column_12").cast(pl.Utf8,  strict=False).alias("LOEUF"),
                  pl.col("column_13").cast(pl.Int64, strict=False).alias("bp_overlap"),
              ])
              # Nettoie et cast LOEUF
              .with_columns([
                  pl.col("gene_type").replace([".", "NA"], None),
                  pl.col("gene_name").replace([".", "NA"], None),
                  pl.col("transcript").replace([".", "NA"], None),
                  pl.col("LOEUF").replace([".", "NA"], None).cast(pl.Float64, strict=False),
              ])
              # gene_id depuis attrs ; crée t_Stop puis supprime attrs et t_End
              .with_columns(
                  pl.when(pl.col("attrs").is_not_null())
                    .then(
                        pl.col("attrs")
                          .str.replace_all('""', '"')
                          .str.strip_chars('"')
                          .str.extract(r'gene_id\s+"([^"]+)\.\d+"')
                    )
                    .otherwise(None)
                    .alias("GeneID")
              )
              .with_columns(pl.col("t_End").alias("t_Stop"))
              .drop("attrs", "t_End")
        )

        # Schéma (ordre) de colonnes côté gène pour reconstituer l’ordre plus tard
        gene_schema_cols = gene_lf.collect_schema().names()

        # CNV x prob_regions (fraction max par CNV)
        preg_lf = (
            pl.scan_csv(preg_inter_bed, separator="\t", has_header=False)
              .select([
                  pl.col("column_1").alias("Chr"),
                  pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                  pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
                  pl.col("column_7").cast(pl.Int64, strict=False).alias("bp_overlap"),
              ])
              .with_columns(
                  (pl.col("bp_overlap") / (pl.col("Stop") - pl.col("Start") + 1.0))
                  .alias("cnv_problematic_region_overlap")
              )
              .group_by(["Chr", "Start", "Stop"])
              .agg(pl.max("cnv_problematic_region_overlap").alias("cnv_problematic_region_overlap"))
        )

        # Joins (lazy → collect)
        joined_lf = (
            cnvLF.join(gene_lf, on=["Chr", "Start", "Stop"], how="left")
                 .join(preg_lf, on=["Chr", "Start", "Stop"], how="left")
        )
        out_df = joined_lf.collect()

        # Déduplication globale (équivalent | sort | uniq)
        out_df = out_df.sort(by=out_df.columns).unique(maintain_order=True)

        # Ordre final imposé
        out_df = _finalize_order(out_df, cnv_extras, gene_schema_cols)

        # Écriture
        out_df.write_csv(args.out, separator="\t")

    finally:
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()
