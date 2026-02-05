#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LEfSe runner (overall + pairwise) using biobakery/LEfSe Python CLI.

Input files (TSV):
- feature_table: rows=ASV/OTU (features), cols=samples, values=counts (raw preferred)
- taxonomy:      columns include the taxonomic ranks (e.g., Kingdom..Species),
                 must contain a column named 'Feature' (or 'ASV','OTU','ID') mapping to feature_table row IDs
- metadata:      rows=samples, must contain the class column (e.g., Group1)

Example:
python run_lefse.py \
  --feature_table GutBiomeTech/MicrobiomeAnalyst/upload/feature_table_for_MA.txt \
  --taxonomy      GutBiomeTech/MicrobiomeAnalyst/upload/taxonomy_for_MA.txt \
  --metadata      GutBiomeTech/MicrobiomeAnalyst/upload/metadata_for_MA.txt \
  --class_col Group1 --rank Genus \
  --outdir GutBiomeTech/MicrobiomeAnalyst/results_LEfSe_py \
  --lda_cut 2.0 --p_cut 0.10
"""

import argparse
import os
import sys
import shutil
import tempfile
from pathlib import Path
import itertools
import pandas as pd
import numpy as np
import subprocess

# ------------------------- helpers -------------------------

def find_lefse_scripts():
    """Return absolute paths to LEfSe scripts or None if missing."""
    names = ["format_input.py", "run_lefse.py", "plot_res.py"]
    paths = []
    for n in names:
        p = shutil.which(n)
        paths.append(p)
    return dict(zip(names, paths))

def read_table_guess(path):
    return pd.read_csv(path, sep="\t", header=0, dtype=str).fillna("")

def coerce_numeric(df):
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.fillna(0.0)
    return out

def detect_feature_id_col(tax):
    candidates = ["Feature","ASV","OTU","ID","feature","asv","otu","id"]
    for c in candidates:
        if c in tax.columns:
            return c
    # fallback: first column
    return tax.columns[0]

def collapse_to_rank(feat_df, tax_df, id_col, rank):
    # feat_df: index=feature id, columns=samples (numeric counts)
    # tax_df : row per feature id, columns include rank
    if rank not in tax_df.columns:
        raise ValueError(f"Rank '{rank}' not found in taxonomy columns: {list(tax_df.columns)}")
    merged = tax_df[[id_col, rank]].copy()
    merged[rank] = merged[rank].replace("", np.nan)
    merged[rank] = merged[rank].fillna(f"Unclassified_{rank}")
    merged = merged.set_index(id_col)
    # align taxonomy to features
    common = merged.index.intersection(feat_df.index)
    if len(common) == 0:
        raise ValueError("No overlap between taxonomy feature IDs and feature table row IDs.")
    feat_aln = feat_df.loc[common]
    rk_labels = merged.loc[common, rank]
    # groupby rank (sum counts within rank)
    feat_aln["_rk"] = rk_labels.values
    collapsed = feat_aln.groupby("_rk").sum(numeric_only=True)
    collapsed.index.name = "Taxon"
    return collapsed

def make_lefse_matrix(collapsed_counts, class_ser):
    """
    Build LEfSe input (row-based class line at the top).
    collapsed_counts: rows=taxa, cols=samples (numeric)
    class_ser: index=samples, values=class labels
    """
    # intersect samples & order
    samples = [s for s in collapsed_counts.columns if s in class_ser.index]
    if len(samples) < 4:
        raise ValueError("Too few samples after matching counts and metadata; need >= 4.")
    X = collapsed_counts.loc[:, samples]
    y = class_ser.loc[samples].astype(str)

    # LEfSe expects a table like:
    # class   samp1 samp2 ...
    # <class>  C    D    ...
    # (no subclass/subject lines when not used)
    # feature rows follow
    top = pd.DataFrame([["class"] + samples, [""] + list(y)], dtype=str)
    # data
    data = X.copy()
    data.insert(0, "Feature", data.index)
    mat = pd.concat([top, data.reset_index(drop=True)], axis=0, ignore_index=True)
    return mat

def run_cmd(cmd, **kwargs):
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
    if r.returncode != 0:
        sys.stderr.write(f"[CMD ERROR] {' '.join(cmd)}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}\n")
        raise RuntimeError("Command failed.")
    return r

def write_index_tsv(outdir):
    out = []
    for root, _, files in os.walk(outdir):
        for f in files:
            if f.lower().endswith((".png",".pdf",".tsv",".res",".in")):
                p = os.path.join(root, f)
                rel = os.path.relpath(p, outdir)
                out.append((rel, p))
    if out:
        idx = pd.DataFrame(out, columns=["relative_path","abs_path"])
        idx.to_csv(os.path.join(outdir, "_index.tsv"), sep="\t", index=False)

# ------------------------- main -------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--feature_table", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--class_col", default="Group1")
    ap.add_argument("--rank", default="Genus")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--p_cut", type=float, default=0.10)
    ap.add_argument("--lda_cut", type=float, default=2.0)  # log10
    ap.add_argument("--min_per_group", type=int, default=3)
    ap.add_argument("--relative", action="store_true", help="Normalize to relative abundance per sample.")
    args = ap.parse_args()

    scripts = find_lefse_scripts()
    missing = [k for k,v in scripts.items() if v is None]
    if missing:
        print(f"[ERROR] Missing LEfSe scripts in PATH: {missing}\n"
              f"Try: conda install -c bioconda -c conda-forge lefse", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "overall").mkdir(exist_ok=True, parents=True)
    (outdir / "pairwise").mkdir(exist_ok=True, parents=True)

    # ---------- load inputs ----------
    ft = read_table_guess(args.feature_table)
    tx = read_table_guess(args.taxonomy)
    md = read_table_guess(args.metadata)

    # feature table: first column may be feature ID if header says so (Feature/ASV/OTU..)
    # detect feature ID column
    # If the first column contains non-numeric, assume it's row IDs.
    ft_num = ft.copy()
    if not np.issubdtype(ft.iloc[:,1].dtype, np.number):
        # try to coerce all but first to numeric
        idx_col = ft.columns[0]
        ft_num = ft.set_index(idx_col)
        ft_num = coerce_numeric(ft_num)
    else:
        # assume first col is a sample, so we can't identify feature IDs -> raise
        print("[WARN] feature_table appears to have samples starting from first column. "
              "Expecting first column to be feature IDs; attempting fallback by making a synthetic ID index.",
              file=sys.stderr)
        ft_num = ft.copy()
        ft_num.index = [f"f{i+1}" for i in range(len(ft))]
        ft_num = coerce_numeric(ft_num)

    id_col = detect_feature_id_col(tx)
    tx[id_col] = tx[id_col].astype(str)
    # make feature IDs string, align
    ft_num.index = ft_num.index.astype(str)

    # metadata: need class column, row index = Sample IDs
    if args.class_col not in md.columns:
        raise ValueError(f"class_col '{args.class_col}' not found in metadata.")
    # try using first column as Sample if present
    if md.columns[0] not in [args.class_col]:
        md = md.set_index(md.columns[0])
    md.index = md.index.astype(str)
    y_all = md[args.class_col].astype(str)

    # ---------- collapse to rank ----------
    collapsed = collapse_to_rank(ft_num, tx, id_col=id_col, rank=args.rank)

    # optional relative normalization
    if args.relative:
        collapsed = collapsed.div(collapsed.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    # save the collapsed matrix
    collapsed_path = outdir / f"collapsed_{args.rank}_taxaXsample_counts.tsv"
    collapsed.assign(Taxon=collapsed.index).to_csv(collapsed_path, sep="\t", index=False)

    # ------------- OVERALL -------------
    overall_dir = outdir / "overall"
    overall_dir.mkdir(exist_ok=True)
    mat_overall = make_lefse_matrix(collapsed, y_all)
    in_txt = overall_dir / "lefse_input.tsv"
    mat_overall.to_csv(in_txt, sep="\t", header=False, index=False)

    with tempfile.TemporaryDirectory() as tmpd:
        formatted_in = overall_dir / "lefse_formatted.in"
        res_path     = overall_dir / "lefse_results.res"
        # format input: -c 1 (class row 1), no subclass/subject
        run_cmd([scripts["format_input.py"], str(in_txt), str(formatted_in),
                 "-c", "1", "-s", "0", "-u", "0"])
        # run lefse
        run_cmd([scripts["run_lefse.py"], str(formatted_in), str(res_path),
                 "-l", str(args.lda_cut),
                 "--kw_cutoff", str(args.p_cut),
                 "--wilcoxon_cutoff", str(args.p_cut),
                 "--fdr"])
        # plot (bar chart)
        run_cmd([scripts["plot_res.py"], str(res_path),
                 str(overall_dir / "lefse_bar.png"),
                 "--format", "png", "--dpi", "200"])

    # also export table as TSV (feature, class, lda)
    # biobakery's res file is tab-delimited (feature, class, lda, p, ...). Dump it.
    try:
        res_txt = pd.read_csv(res_path, sep="\t", header=None)
        res_txt.to_csv(overall_dir / "lefse_results.tsv", sep="\t", index=False,
                       header=["feature","class","lda","pvalue","something"])
    except Exception:
        pass

    # ------------- PAIRWISE -------------
    classes = sorted(y_all.dropna().unique())
    pairs = list(itertools.combinations(classes, 2))
    for (g1, g2) in pairs:
        sub_dir = outdir / "pairwise" / f"{g1}_vs_{g2}"
        sub_dir.mkdir(parents=True, exist_ok=True)
        # filter samples
        sel_samps = y_all.index[y_all.isin([g1,g2])]
        sub_y = y_all.loc[sel_samps]
        if min((sub_y==g1).sum(), (sub_y==g2).sum()) < args.min_per_group:
            # too few per group
            with open(sub_dir / "SKIPPED.txt","w") as f:
                f.write(f"Skipped {g1} vs {g2}: need >= {args.min_per_group} per group.\n")
            continue
        sub_mat = make_lefse_matrix(collapsed.loc[:, collapsed.columns.isin(sel_samps)], sub_y)
        in_txt = sub_dir / "lefse_input.tsv"
        sub_mat.to_csv(in_txt, sep="\t", header=False, index=False)
        formatted_in = sub_dir / "lefse_formatted.in"
        res_path     = sub_dir / "lefse_results.res"

        run_cmd([scripts["format_input.py"], str(in_txt), str(formatted_in),
                 "-c", "1", "-s", "0", "-u", "0"])
        run_cmd([scripts["run_lefse.py"], str(formatted_in), str(res_path),
                 "-l", str(args.lda_cut),
                 "--kw_cutoff", str(args.p_cut),
                 "--wilcoxon_cutoff", str(args.p_cut),
                 "--fdr"])
        run_cmd([scripts["plot_res.py"], str(res_path),
                 str(sub_dir / "lefse_bar.png"),
                 "--format", "png", "--dpi", "200"])

        # TSV export (best-effort)
        try:
            res_txt = pd.read_csv(res_path, sep="\t", header=None)
            res_txt.to_csv(sub_dir / "lefse_results.tsv", sep="\t", index=False,
                           header=["feature","class","lda","pvalue","something"])
        except Exception:
            pass

    write_index_tsv(outdir)
    print(f"== Done ==\nOutput: {outdir}\n- overall/* and pairwise/* created.")

if __name__ == "__main__":
    main()

