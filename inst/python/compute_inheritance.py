#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CNV inheritance annotation (child vs parents) using reciprocal overlap.

Inputs
------
- --cnv_geneAnnot : Gene-annotated CNV table (TSV / TSV.GZ), e.g. columns
  SampleID, Chr, Start, End or Stop, Type, gene_name.
- --pedigree      : Pedigree table (TSV) containing ChildID / FatherID / MotherID
  (column names auto-detected with sensible fallbacks).
- --overlap       : Single reciprocal overlap fraction threshold (child and parent),
  e.g. 0.5 means >=50% overlap on both the child CNV and the parent CNV.

Output
------
- --output        : TSV containing CHILD rows only, with:
    - transmitted_cnv  (bool): True=inherited, False=de novo
    - transmitted_gene (bool / "intergenic" / NaN)

Notes
-----
- Overlap is computed with bioframe.overlap (0-based half-open intervals).
- CNV type matching is enforced (DEL with DEL, DUP with DUP).
- End is auto-harmonized from Stop/STOP if needed.
- If a parent belongs to multiple TrioKeys (multiple probands), a many-to-many
  merge is allowed so each link is evaluated.
"""

from __future__ import annotations
import argparse
import sys
import pandas as pd
import numpy as np

# Dependency: bioframe
try:
    import bioframe as bf
except Exception as e:
    print("Error: this script requires 'bioframe' (pip install bioframe).", file=sys.stderr)
    raise


# ============================ Utilities ============================

def read_table_auto(path: str) -> pd.DataFrame:
    """
    Read a tab-delimited file (TSV/TSV.GZ). Uses pandas' default C engine.
    """
    return pd.read_table(path, sep="\t", dtype=str, low_memory=False)


def _find_col(df: pd.DataFrame, candidates) -> str | None:
    """
    Return the first matching column name among `candidates` using case-insensitive
    equality or startswith() matching. Returns None if nothing matches.
    """
    cols = list(df.columns)
    lower_map = {c: c.lower() for c in cols}
    for cand in candidates:
        for c in cols:
            lc = lower_map[c]
            if lc == cand or lc.startswith(cand):
                return c
    return None


def load_and_normalize_pedigree(path_or_df, sep: str = r"\t") -> pd.DataFrame:
    """
    Normalize pedigree to columns: ChildID, FatherID, MotherID (dtype 'string').
    Accepts a DataFrame or a path.
    """
    if isinstance(path_or_df, pd.DataFrame):
        ped_raw = path_or_df.copy()
    else:
        ped_raw = pd.read_csv(path_or_df, sep=sep, dtype=str)

    child_col  = _find_col(ped_raw, ["sampleid", "child", "proband", "subject"])
    father_col = _find_col(ped_raw, ["father", "biofather", "fatherid", "biofatherid", "dad"])
    mother_col = _find_col(ped_raw, ["mother", "biomother", "motherid", "biomotherid", "mom"])

    if child_col  is None and "SampleID" in ped_raw.columns: child_col  = "SampleID"
    if father_col is None and "FatherID"  in ped_raw.columns: father_col = "FatherID"
    if mother_col is None and "MotherID"  in ped_raw.columns: mother_col = "MotherID"

    rename_map = {}
    if child_col:  rename_map[child_col]  = "ChildID"
    if father_col: rename_map[father_col] = "FatherID"
    if mother_col: rename_map[mother_col] = "MotherID"

    ped = ped_raw.rename(columns=rename_map)
    for c in ["ChildID", "FatherID", "MotherID"]:
        if c not in ped.columns:
            ped[c] = pd.NA

    ped = ped[["ChildID", "FatherID", "MotherID"]].astype("string")
    for c in ["ChildID", "FatherID", "MotherID"]:
        ped[c] = ped[c].astype("string").str.strip()

    return ped


def keep_complete_trios(ped: pd.DataFrame) -> pd.DataFrame:
    """
    Keep complete trios only and create a unique TrioKey.
    """
    m = (~ped["ChildID"].isna() & (ped["ChildID"] != "") &
         ~ped["FatherID"].isna() & (ped["FatherID"] != "") &
         ~ped["MotherID"].isna() & (ped["MotherID"] != ""))
    out = ped.loc[m].copy()
    out["TrioKey"] = out["ChildID"] + "_" + out["FatherID"] + "_" + out["MotherID"]
    out = out.drop_duplicates(subset=["TrioKey"]).reset_index(drop=True)
    return out[["ChildID", "FatherID", "MotherID", "TrioKey"]]


def canonicalize_cnv_columns(cnv: pd.DataFrame) -> pd.DataFrame:
    """
    Detect and rename CNV columns to canonical names:
      SampleID, Chr, Start, End, Type, (gene_name optional)

    This makes the script robust to common header variants like:
    - stop/STOP/END -> End
    - chr/CHROM/chromosome -> Chr
    - svtype/type -> Type
    - symbol/gene -> gene_name
    """
    col_map = {}

    sample = _find_col(cnv, ["sampleid", "sample", "iid", "id"])
    chrom  = _find_col(cnv, ["chr", "chrom", "chromosome", "contig"])
    start  = _find_col(cnv, ["start", "pos", "begin", "position"])
    end    = _find_col(cnv, ["end", "stop"])
    typ    = _find_col(cnv, ["type", "svtype"])
    gene   = _find_col(cnv, ["gene_name", "gene", "symbol"])

    if sample is None and "SampleID" in cnv.columns: sample = "SampleID"
    if chrom  is None and "Chr"      in cnv.columns: chrom  = "Chr"
    if start  is None and "Start"    in cnv.columns: start  = "Start"
    if end    is None and "End"      in cnv.columns: end    = "End"
    if typ    is None and "Type"     in cnv.columns: typ    = "Type"
    # gene is optional

    if sample: col_map[sample] = "SampleID"
    if chrom:  col_map[chrom]  = "Chr"
    if start:  col_map[start]  = "Start"
    if end:    col_map[end]    = "End"
    if typ:    col_map[typ]    = "Type"
    if gene:   col_map[gene]   = "gene_name"

    cnv2 = cnv.rename(columns=col_map)

    required = ["SampleID", "Chr", "Start", "End", "Type"]
    missing = [c for c in required if c not in cnv2.columns]
    if missing:
        raise ValueError(f"Missing columns in CNV (after normalization): {missing}")

    return cnv2


# ============================ Pipeline ============================

def run_pipeline(
    cnv_geneannot_path: str,
    pedigree_path: str,
    output_path: str,
    overlap_frac: float,
):
    # --- Pedigree: load, normalize, keep complete trios
    ped_src  = read_table_auto(pedigree_path)
    ped_norm = load_and_normalize_pedigree(ped_src, sep=r"\t")

    # Guard: check that pedigree columns were recognized
    all_na = ped_norm[["ChildID", "FatherID", "MotherID"]].isna().all().all()
    if all_na:
        raise ValueError(
            "Pedigree file loaded but no ChildID/FatherID/MotherID columns could be detected. "
            "Expected column names (case-insensitive): SampleID/Child/Proband for child, "
            "Father/FatherID/Dad for father, Mother/MotherID/Mom for mother."
        )

    trios = keep_complete_trios(ped_norm).replace({'': pd.NA})

    # Guard: no complete trios found
    if trios.empty:
        raise ValueError(
            f"No complete trios found in pedigree (rows with non-empty ChildID + FatherID + MotherID). "
            f"Pedigree has {len(ped_norm)} rows. Check that all three parent columns are filled."
        )
    print(f"INFO: {len(trios)} complete trio(s) found.", file=sys.stderr)

    # Long format: (TrioKey, SampleID, family_status)
    role_map = {"ChildID": "child", "FatherID": "father", "MotherID": "mother"}
    trios_reformated = (
        trios.melt(
            id_vars=["TrioKey"],
            value_vars=["ChildID", "FatherID", "MotherID"],
            var_name="family_status",
            value_name="SampleID",
        )
        .dropna(subset=["SampleID"])
        .assign(family_status=lambda d: d["family_status"].map(role_map))
        [["TrioKey", "SampleID", "family_status"]]
        .drop_duplicates(subset=["TrioKey", "SampleID", "family_status"])
    )

    # --- CNVs: load and harmonize expected columns
    cnv_raw = read_table_auto(cnv_geneannot_path)
    cnv     = canonicalize_cnv_columns(cnv_raw)

    # Guard: CNV file is empty
    if cnv.empty:
        raise ValueError("CNV file is empty (no rows after loading).")
    print(f"INFO: {len(cnv)} CNV row(s) loaded.", file=sys.stderr)

    # Join CNVs with trio roles (allow many-to-many to support parents shared across trios)
    merged = pd.merge(
        cnv,
        trios_reformated,
        on="SampleID",
        how="inner"
    )

    # Guard: no overlap between CNV SampleIDs and pedigree SampleIDs
    if merged.empty:
        cnv_ids  = set(cnv["SampleID"].dropna().unique())
        ped_ids  = set(trios_reformated["SampleID"].dropna().unique())
        common   = cnv_ids & ped_ids
        raise ValueError(
            f"No CNV rows matched any sample in the pedigree (inner join returned 0 rows). "
            f"CNV file has {len(cnv_ids)} unique SampleIDs, pedigree has {len(ped_ids)}, "
            f"common IDs: {len(common)}. Check that SampleID values match between both files."
        )

    # ===== CNV-level inheritance (segment-based) =====
    VALID_TYPES = {"DEL", "DUP"}

    df = merged.copy()

    # Normalize Type and family roles
    df["Type"] = (
        df["Type"].astype(str).str.strip().str.upper()
          .replace({"DELETION": "DEL", "DUPLICATION": "DUP"})
    )
    df["family_status"] = df["family_status"].astype(str).str.strip().str.lower()

    # Coordinates as numeric
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"]   = pd.to_numeric(df["End"],   errors="coerce")

    # Basic filtering
    df = df[df["Type"].isin(VALID_TYPES)].copy()
    df = df.dropna(subset=["TrioKey", "Chr", "Start", "End", "Type", "family_status"])
    df = df.query("End > Start").copy()

    # Child CNV ID (used to propagate inheritance to all rows)
    df["cnv_id"] = (
        df["TrioKey"].astype(str) + "_" +
        df["Chr"].astype(str) + "_" +
        df["Start"].astype("Int64").astype(str) + "_" +
        df["End"].astype("Int64").astype(str) + "_" +
        df["Type"]
    )

    is_child  = df["family_status"].eq("child")
    is_parent = df["family_status"].isin(["father", "mother"])

    child_uni = (
        df.loc[is_child, ["TrioKey", "Chr", "Start", "End", "Type", "cnv_id"]]
          .drop_duplicates()
          .rename(columns={"Chr": "chrom", "Start": "start", "End": "end"})
    )
    parents = (
        df.loc[is_parent, ["TrioKey", "Chr", "Start", "End", "Type"]]
          .drop_duplicates()
          .rename(columns={"Chr": "chrom", "Start": "start", "End": "end"})
    )

    # Reciprocal overlap (same TrioKey & Type) via bioframe
    if not child_uni.empty and not parents.empty:
        ov = bf.overlap(
            child_uni[["chrom", "start", "end", "TrioKey", "Type", "cnv_id"]],
            parents  [["chrom", "start", "end", "TrioKey", "Type"]],
            on=["TrioKey", "Type"],
            how="inner",
            suffixes=("", "_p"),
        )
        if not ov.empty:
            # overlap_len = max(0, min(end, end_p) - max(start, start_p))
            ov_len    = np.clip(np.minimum(ov["end"], ov["end_p"]) - np.maximum(ov["start"], ov["start_p"]), a_min=0, a_max=None)
            child_len = (ov["end"]   - ov["start"]).clip(lower=1)
            parent_len= (ov["end_p"] - ov["start_p"]).clip(lower=1)

            thr = float(overlap_frac)
            ok_child  = (ov_len / child_len)  >= thr
            ok_parent = (ov_len / parent_len) >= thr

            inherited_ids = set(ov.loc[ok_child & ok_parent, "cnv_id"])
        else:
            inherited_ids = set()
    else:
        inherited_ids = set()

    child_status = pd.DataFrame({
        "cnv_id": child_uni["cnv_id"],
        "transmitted_cnv": child_uni["cnv_id"].isin(inherited_ids)  # True=inherited, False=de novo
    })

    df = df.merge(child_status, on="cnv_id", how="left")

    # ===== Gene-level (per-row): True / False / "intergenic" =====
    gdf = df.copy()

    gdf["Type"]          = gdf["Type"].astype(str).str.strip().str.upper()
    gdf["family_status"] = gdf["family_status"].astype(str).str.strip().str.lower()

    if "gene_name" not in gdf.columns:
        gdf["transmitted_gene"] = pd.NA
    else:
        gdf["gene_name_norm"] = (
            gdf["gene_name"].astype("string").str.strip()
               .replace({'': pd.NA, 'nan': pd.NA, 'NaN': pd.NA, 'None': pd.NA})
        )

        is_child_g  = gdf["family_status"].eq("child")
        is_parent_g = gdf["family_status"].isin(["father", "mother"])
        has_gene    = gdf["gene_name_norm"].notna()

        parent_gene_keys = (
            gdf.loc[is_parent_g & has_gene, ["TrioKey", "gene_name_norm", "Type"]]
              .drop_duplicates()
              .assign(_parent_has=True)
        )
        child_gene_keys = (
            gdf.loc[is_child_g & has_gene, ["TrioKey", "gene_name_norm", "Type"]]
              .drop_duplicates()
        )

        tg = child_gene_keys.merge(
                parent_gene_keys,
                on=["TrioKey", "gene_name_norm", "Type"],
                how="left"
             )
        # Boolean for children with a gene in the row: True=inherited (present in a parent), False=de novo
        tg["transmitted_gene"] = tg["_parent_has"].fillna(False)
        tg = tg.drop(columns=["_parent_has"])

        gdf = gdf.merge(tg, on=["TrioKey", "gene_name_norm", "Type"], how="left")
        # Intergenic rows for children (no gene on the row)
        gdf.loc[is_child_g & (~has_gene), "transmitted_gene"] = "intergenic"
        # Hide for parents
        gdf.loc[~is_child_g, "transmitted_gene"] = np.nan
        gdf = gdf.drop(columns=["gene_name_norm"])

    # ===== Output: CHILD rows only =====
    child_df = gdf.loc[gdf["family_status"].eq("child")].copy()

    # Guard: no child rows in output
    if child_df.empty:
        raise ValueError(
            "No child rows found after processing. "
            "This can happen if: (1) no CNVs belong to child samples, "
            "(2) all child CNVs were filtered out (invalid Type, missing coordinates), or "
            "(3) SampleIDs in the CNV file do not match ChildIDs in the pedigree."
        )

    n_inherited = child_df["transmitted_cnv"].eq(True).sum()
    n_denovo    = child_df["transmitted_cnv"].eq(False).sum()
    print(
        f"INFO: {len(child_df)} child row(s) — "
        f"{n_inherited} inherited CNV(s), {n_denovo} de novo CNV(s).",
        file=sys.stderr
    )
    if n_inherited == 0:
        print(
            "WARN: No inherited CNVs found. All child CNVs are marked de novo. "
            "Check overlap threshold (--overlap) and that parent CNVs are present in the input file.",
            file=sys.stderr
        )

    child_df.to_csv(output_path, sep="\t", index=False)
    print(f"OK: {len(child_df)} CHILD rows written -> {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate CNV inheritance (child vs parents) using reciprocal overlap; 'Type' column is assumed."
    )
    parser.add_argument("--cnv_geneAnnot", required=True,
                        help="Gene-annotated CNV file (TSV/TSV.GZ). Expected columns (or variants auto-detected): SampleID, Chr, Start, End/Stop, Type, [gene_name].")
    parser.add_argument("--pedigree", required=True,
                        help="Pedigree file (TSV) with ChildID/FatherID/MotherID (column names auto-detected).")
    parser.add_argument("--output", required=True,
                        help="Output TSV (CHILD rows only).")
    parser.add_argument("--overlap", type=float, default=0.5,
                        help="Single reciprocal overlap threshold for child and parent. Default: 0.5")

    args = parser.parse_args()

    try:
        run_pipeline(
            cnv_geneannot_path=args.cnv_geneAnnot,
            pedigree_path=args.pedigree,
            output_path=args.output,
            overlap_frac=args.overlap,
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
