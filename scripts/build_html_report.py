#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, glob, shutil
from pathlib import Path
from collections import defaultdict

# ---------- Styles ----------
CSS = """
:root{--bg:#0b1020;--card:#111633;--ink:#e6ecff;--muted:#aab3d7;--accent:#8aa1ff}
*{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--ink);font-family:ui-sans-serif,system-ui,Segoe UI,Roboto,Helvetica,Arial}
a{color:var(--accent);text-decoration:none}
.wrapper{display:grid;grid-template-columns:260px 1fr;min-height:100vh}
nav{background:#0d1330;border-right:1px solid #1b2457;padding:20px;position:sticky;top:0;height:100vh;overflow:auto}
main{padding:28px}
h1{margin:0 0 16px 0;font-size:22px}
h2{font-size:20px;margin:24px 0 12px}
.section{margin-bottom:34px}
.card{background:var(--card);border:1px solid #1b2457;border-radius:16px;padding:16px;margin:12px 0}
.list a{display:block;padding:10px 12px;border:1px solid #27306a;background:#0e1433;border-radius:10px;margin:8px 0}
.breadcrumbs{font-size:12px;color:var(--muted);margin:4px 0 16px}
.grid{display:grid;grid-template-columns:repeat(auto-fill, minmax(280px, 1fr));gap:12px}
.thumb{background:#0e1433;border:1px solid #27306a;border-radius:12px;padding:10px}
.thumb .name{font-size:12px;color:var(--muted);margin-top:6px;word-break:break-all}
.thumb .src{font-size:12px;margin-top:6px}
.badge{display:inline-block;padding:2px 8px;border-radius:999px;background:#1a2259;border:1px solid #2a3580;color:#c8d3ff;font-size:12px}
.count{font-size:12px;color:#aab3d7;margin-left:8px}
.note{color:var(--muted);font-size:13px}
hr{border:0;border-top:1px solid #1b2457;margin:24px 0}
footer{color:#8590c8;font-size:12px;margin-top:24px}
kbd{background:#1a2259;color:#c8d3ff;border-radius:6px;padding:2px 6px;border:1px solid #2a3580}
"""

LAYOUT_HEAD = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>{title}</title><meta name="viewport" content="width=device-width, initial-scale=1">
<style>{css}</style></head>
<body><div class="wrapper"><nav>
<h1>🧫 Microbiome Report</h1>
<ol>
  <li><a href="{root}index.html">Overview</a></li>
  <li><a href="{root}composition/index.html">Composition</a></li>
  <li><a href="{root}alpha/index.html">Alpha diversity</a></li>
  <li><a href="{root}beta/index.html">Beta diversity</a></li>
  <li><a href="{root}lefse/index.html">LEfSe</a></li>
  <li><a href="{root}picrust/index.html">PICRUSt2 Pathway</a></li>
</ol>
</nav><main>
<div class="breadcrumbs">{breadcrumbs}</div>
"""

LAYOUT_FOOT = """
<footer>Self-contained report • Generated locally</footer>
</main></div></body></html>
"""

IMG_EXT = (".png", ".jpg", ".jpeg", ".gif", ".svg")
PDF_EXT = (".pdf",)

# MA_OUTDIR 표준 경로
TAX_LEVELS = ["phylum","class","order","family","genus","species"]
ALPHA_TSV = "01_Alpha/alpha_values.tsv"
BETA_MAIN_PAT = "02_Beta/perma*bray_all.tsv"
BETA_PAIR = "02_Beta/pairwise_permanova_bray.tsv"
BETA_PCOA_SRC = "02_Beta/pcoa_bray_genus_source.tsv"
COMP_PROPS_DIR = "03_Composition/props"

# ---------- Helpers ----------
def list_files(root, pats=("**/*.png","**/*.pdf","**/*.svg","**/*.jpg","**/*.jpeg")):
    files = []
    for p in pats:
        files += glob.glob(os.path.join(root, p), recursive=True)
    return sorted(set(files))

def materialize(src_files, dest_dir, mode="copy"):
    dest_dir = Path(dest_dir); dest_dir.mkdir(parents=True, exist_ok=True)
    mapped = []
    for src in src_files:
        sp = Path(src)
        if not sp.exists():
            continue
        dst = dest_dir / sp.name
        i = 1
        while dst.exists():
            dst = dest_dir / f"{sp.stem}__{i}{sp.suffix}"
            i += 1
        if mode == "hardlink":
            try: os.link(sp.resolve(), dst)
            except OSError: shutil.copy2(sp, dst)
        else:
            shutil.copy2(sp, dst)
        mapped.append(dst.resolve())
    return mapped

def rels_from(page_dir: Path, abs_paths):
    page_dir = Path(page_dir).resolve()
    return [os.path.relpath(Path(p).resolve(), start=page_dir) for p in abs_paths]

def write_html(out_path: Path, title: str, root_href: str, breadcrumbs_html: str, inner_html: str):
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    html = LAYOUT_HEAD.format(title=title, css=CSS, root=root_href, breadcrumbs=breadcrumbs_html) + inner_html + LAYOUT_FOOT
    out_path.write_text(html, encoding="utf-8")

def make_list_page(title, items):
    parts = [f'<section class="section"><h2>{title}</h2><div class="card list">']
    for label, href, note in items:
        tail = f' <span class="count">({note})</span>' if note not in (None,"") else ""
        parts.append(f'<a href="{href}">{label}{tail}</a>')
    parts.append('</div></section>')
    return "\n".join(parts)

def gallery(title, page_dir: Path, asset_hrefs, with_source_links=True):
    if not asset_hrefs:
        return f'<section class="section"><h2>{title}<span class="count">(0 files)</span></h2><div class="note">No files.</div></section>'
    parts = [f'<section class="section"><h2>{title}<span class="count">({len(asset_hrefs)} files)</span></h2>']
    imgs = [h for h in asset_hrefs if h.lower().endswith(IMG_EXT)]
    pdfs = [h for h in asset_hrefs if h.lower().endswith(PDF_EXT)]

    if imgs:
        parts += ['<div class="card"><div class="badge">Images</div><div class="grid">']
        for h in imgs:
            src_html = ""
            if with_source_links:
                abs_img = (page_dir / h).resolve()
                src_abs = abs_img.with_name(re.sub(r'\.(png|jpg|jpeg|gif|svg)$', '_source.tsv', abs_img.name, flags=re.I))
                if src_abs.exists():
                    src_rel = os.path.relpath(src_abs, start=page_dir)
                    src_html = f'<div class="src">Data: <a href="{src_rel}" target="_blank">TSV</a></div>'
            parts += [
                f'<div class="thumb">'
                f'<a href="{h}" target="_blank"><img src="{h}" style="width:100%;border-radius:8px"/></a>'
                f'<div class="name">{h}</div>{src_html}'
                f'</div>'
            ]
        parts += ['</div></div>']

    if pdfs:
        parts += ['<div class="card"><div class="badge">PDFs</div>']
        for h in pdfs:
            parts += [f'<div class="thumb"><div class="name">{h}</div>'
                      f'<object data="{h}" type="application/pdf" width="100%" height="600px">'
                      f'<p><a href="{h}" target="_blank">Open</a></p></object></div>']
        parts += ['</div>']

    parts += ['</section>']
    return "\n".join(parts)

def data_card(title, page_dir: Path, abs_files):
    if not abs_files:
        return ""
    hrefs = rels_from(page_dir, abs_files)
    items = "".join([f'<li><a href="{h}" target="_blank">{Path(h).name}</a></li>' for h in hrefs])
    return f'<section class="section"><h2>{title}</h2><div class="card"><ul>{items}</ul></div></section>'

def split_by_tax(files):
    buckets = {lvl: [] for lvl in TAX_LEVELS}
    for f in files:
        lname = f.lower(); matched = False
        for lvl in TAX_LEVELS:
            if lvl in lname:
                buckets[lvl].append(f); matched = True; break
        if not matched:
            buckets["genus"].append(f)
    return buckets

# ---------- LEfSe 보완 스캐너 ----------
def scan_lefse_images(lefse_root: Path):
    out = {
        "overall_bar": [],
        "overall_dot": [],
        "pairwise": defaultdict(lambda: {"bar": [], "dot": []}),
        "other": [],
    }
    if not lefse_root or not lefse_root.exists():
        return out

    std_overall_bar = sorted(glob.glob(str(lefse_root / "overall" / "bar" / "*.png")))
    std_overall_dot = sorted(glob.glob(str(lefse_root / "overall" / "dot" / "*.png")))
    out["overall_bar"] += std_overall_bar
    out["overall_dot"] += std_overall_dot

    for comp_dir in sorted(glob.glob(str(lefse_root / "pairwise" / "*"))):
        comp = Path(comp_dir).name
        bars = sorted(glob.glob(str(Path(comp_dir) / "bar" / "*.png")))
        dots = sorted(glob.glob(str(Path(comp_dir) / "dot" / "*.png")))
        if bars or dots:
            out["pairwise"][comp]["bar"] += bars
            out["pairwise"][comp]["dot"] += dots

    all_pngs = set(glob.glob(str(lefse_root / "**" / "*.png"), recursive=True))
    already = set(out["overall_bar"]) | set(out["overall_dot"])
    for comp, dd in out["pairwise"].items():
        already |= set(dd["bar"]) | set(dd["dot"])
    remain = [p for p in all_pngs if p not in already]

    bar_pat = re.compile(r"(?:^|[_\-./])bar(?:[_\-./]|$)", re.IGNORECASE)
    dot_pat = re.compile(r"(?:^|[_\-./])dot(?:[_\-./]|$)", re.IGNORECASE)
    comp_pat = re.compile(r"([A-Za-z0-9.+]+)[_\-./]*v[sS][_\-./]*([A-Za-z0-9.+]+)")

    for p in remain:
        lp = p.lower()
        m = comp_pat.search(lp)
        if m:
            comp = f"{m.group(1)}_vs_{m.group(2)}"
            if bar_pat.search(lp):
                out["pairwise"][comp]["bar"].append(p); continue
            if dot_pat.search(lp):
                out["pairwise"][comp]["dot"].append(p); continue
        if "overall" in lp:
            if bar_pat.search(lp):
                out["overall_bar"].append(p); continue
            if dot_pat.search(lp):
                out["overall_dot"].append(p); continue
        if bar_pat.search(lp):
            out["overall_bar"].append(p); continue
        if dot_pat.search(lp):
            out["overall_dot"].append(p); continue
        out["other"].append(p)

    out["overall_bar"] = sorted(set(out["overall_bar"]))
    out["overall_dot"] = sorted(set(out["overall_dot"]))
    for comp in list(out["pairwise"].keys()):
        out["pairwise"][comp]["bar"] = sorted(set(out["pairwise"][comp]["bar"]))
        out["pairwise"][comp]["dot"] = sorted(set(out["pairwise"][comp]["dot"]))
    out["other"] = sorted(set(out["other"]))
    return out

# ---------- PICRUSt2 스캐너 (PDF 제한: errorbar_by_class.pdf만) ----------
def scan_picrust_dir(pic_plot_dir: Path):
    """
    plots_by_class 디렉토리에서 PDF는 오직 errorbar_by_class.pdf만,
    TSV는 소스/Top terms/DAA 3종을 모아 반환.
    """
    out = {"figs": [], "data": []}
    if not pic_plot_dir or not pic_plot_dir.exists():
        return out

    # ★ PDF 제한
    only_pdf = pic_plot_dir / "errorbar_by_class.pdf"
    if only_pdf.exists():
        out["figs"] = [str(only_pdf)]

    # TSV들
    tsv_candidates = [
        "errorbar_by_class_source.tsv",
        "errorbar_by_class_top_terms.tsv",
        "daa_kegg_raw.tsv",
        "daa_kegg_kruskal.tsv",
        "daa_kegg_annotated.tsv",
    ]
    data = []
    for name in tsv_candidates:
        p = pic_plot_dir / name
        if p.exists():
            data.append(str(p))
    out["data"] = sorted(set(data))
    return out

# ---------- Builder ----------
def build_site(ma_outdir, ma_fig_dir, lefse_csv, lefse_dir, pic_plot_dir, outdir, mode="copy"):
    outdir = Path(outdir).resolve(); assets = outdir / "assets"
    ma_outdir = Path(ma_outdir).resolve()

    # 1) 그림 스캔 (MA)
    ma_files = list_files(ma_fig_dir)

    # 2) MA_OUTDIR 데이터 파일 경로
    alpha_tsv = ma_outdir / ALPHA_TSV
    beta_main = sorted(ma_outdir.glob(BETA_MAIN_PAT))
    beta_pair = ma_outdir / BETA_PAIR
    beta_pcoa = ma_outdir / BETA_PCOA_SRC
    comp_props = ma_outdir / COMP_PROPS_DIR

    # 3) Composition
    comp_group_all = [f for f in ma_files if ("samgrpbars" in f.lower()) and ("groupmean" in f.lower()) and f.lower().endswith(".png")]
    comp_sample = [f for f in ma_files if "samgrpbars" not in f.lower()
                   and any(w in f.lower() for w in ["composition","barplot","stack","taxa","abundance"])]
    comp_group_tax = split_by_tax(comp_group_all)
    comp_sample_tax = split_by_tax(comp_sample)

    # 4) Alpha/Beta 이미지
    alpha_imgs = [f for f in ma_files if "alpha" in f.lower()]
    beta_imgs = [f for f in ma_files if any(w in f.lower() for w in ["beta","pcoa","nmds","ordination","bray","jaccard","unifrac","dlrp"])]

    # ---------- Overview ----------
    write_html(outdir/"index.html", "Overview", root_href="", breadcrumbs_html="Home",
               inner_html=make_list_page("Choose a section", [
                   ("Composition", "composition/index.html", None),
                   ("Alpha diversity", "alpha/index.html", None),
                   ("Beta diversity", "beta/index.html", None),
                   ("LEfSe", "lefse/index.html", None),
                   ("PICRUSt2 Pathway", "picrust/index.html", None),
               ]))

    # ---------- Composition ----------
    write_html(outdir/"composition/index.html", "Composition", root_href="../",
               breadcrumbs_html='Home / Composition',
               inner_html=make_list_page("Choose view", [
                   ("Group view (groupmean only)", "group/index.html", sum(len(v) for v in comp_group_tax.values())),
                   ("Sample view", "sample/index.html", sum(len(v) for v in comp_sample_tax.values())),
               ]))

    # Group → level list
    grp_level_items = []
    for lvl in TAX_LEVELS:
        cnt = len(comp_group_tax[lvl])
        if cnt > 0: grp_level_items.append((lvl.capitalize(), f"{lvl}.html", cnt))
    write_html(outdir/"composition/group/index.html", "Composition • Group", root_href="../../",
               breadcrumbs_html='Home / Composition / Group',
               inner_html=make_list_page("Choose taxonomic level", grp_level_items))

    # Group leaf pages
    for lvl in TAX_LEVELS:
        files = comp_group_tax[lvl]
        if not files: continue
        srcs = []
        for f in files:
            cand = re.sub(r'\.png$', '_source.tsv', f, flags=re.I)
            if Path(cand).exists(): srcs.append(cand)
        mapped_all = materialize(files + srcs, assets / f"composition/group/{lvl}", mode)
        mapped_imgs = [p for p in mapped_all if str(p).lower().endswith(".png")]
        page_dir = outdir / "composition" / "group"
        hrefs = rels_from(page_dir, mapped_imgs)
        data_files = []
        if comp_props.is_dir():
            base_names = [
                f"{lvl}_wide_samplexTaxon.tsv",
                f"{lvl}_group_mean.tsv",
                f"{lvl}_group_median.tsv",
                f"{lvl}_long.tsv",
                f"{lvl}_taxa_totals.tsv",
            ]
            for bn in base_names:
                cand = comp_props / bn
                if cand.exists():
                    data_files.append(cand)
                else:
                    cand2 = comp_props / bn.replace(lvl, lvl.capitalize())
                    if cand2.exists(): data_files.append(cand2)
        mapped_data = materialize(data_files, assets / f"composition/group/{lvl}", mode)
        inner = gallery(f"{lvl.capitalize()} (group mean only)", page_dir, hrefs) + \
                data_card("Data sources (props)", page_dir, mapped_data)
        write_html(outdir/f"composition/group/{lvl}.html",
                   f"Composition • Group • {lvl.capitalize()}",
                   root_href="../../../",
                   breadcrumbs_html=f'Home / Composition / Group / {lvl.capitalize()}',
                   inner_html=inner)

    # Sample → level list
    smp_level_items = []
    for lvl in TAX_LEVELS:
        cnt = len(comp_sample_tax[lvl])
        if cnt > 0: smp_level_items.append((lvl.capitalize(), f"{lvl}.html", cnt))
    write_html(outdir/"composition/sample/index.html", "Composition • Sample", root_href="../../",
               breadcrumbs_html='Home / Composition / Sample',
               inner_html=make_list_page("Choose taxonomic level", smp_level_items))

    # Sample leaf pages
    for lvl in TAX_LEVELS:
        files = comp_sample_tax[lvl]
        if not files: continue
        mapped = materialize(files, assets / f"composition/sample/{lvl}", mode)
        page_dir = outdir / "composition" / "sample"
        hrefs = rels_from(page_dir, mapped)
        html = gallery(f"{lvl.capitalize()}", page_dir, hrefs, with_source_links=False)
        write_html(outdir/f"composition/sample/{lvl}.html",
                   f"Composition • Sample • {lvl.capitalize()}",
                   root_href="../../../",
                   breadcrumbs_html=f'Home / Composition / Sample / {lvl.capitalize()}',
                   inner_html=html)

    # ---------- Alpha ----------
    mapped_alpha_imgs = materialize(alpha_imgs, assets / "alpha/imgs", mode)
    page_dir_alpha = outdir / "alpha"
    hrefs_alpha = rels_from(page_dir_alpha, mapped_alpha_imgs)
    alpha_data_mapped = materialize([alpha_tsv], assets / "alpha/data", mode) if alpha_tsv.exists() else []
    write_html(outdir/"alpha/index.html", "Alpha diversity", root_href="../",
               breadcrumbs_html='Home / Alpha diversity',
               inner_html=gallery("Alpha diversity plots", page_dir_alpha, hrefs_alpha, with_source_links=False) +
                          data_card("Data sources", page_dir_alpha, alpha_data_mapped))

    # ---------- Beta ----------
    mapped_beta_imgs = materialize(beta_imgs, assets / "beta/imgs", mode)
    page_dir_beta = outdir / "beta"
    hrefs_beta = rels_from(page_dir_beta, mapped_beta_imgs)
    beta_data = []
    beta_data += beta_main
    if beta_pair.exists(): beta_data.append(beta_pair)
    if beta_pcoa.exists(): beta_data.append(beta_pcoa)
    mapped_beta_data = materialize(beta_data, assets / "beta/data", mode) if beta_data else []
    write_html(outdir/"beta/index.html", "Beta diversity", root_href="../",
               breadcrumbs_html='Home / Beta diversity',
               inner_html=gallery("Beta diversity plots", page_dir_beta, hrefs_beta, with_source_links=False) +
                          data_card("Data sources", page_dir_beta, mapped_beta_data))

    # ---------- LEfSe (이전 스타일) ----------
    page_dir_lefse = outdir / "lefse"; os.makedirs(page_dir_lefse, exist_ok=True)
    lefse_root = None
    if lefse_dir and Path(lefse_dir).exists():
        lefse_root = Path(lefse_dir).resolve()
    elif lefse_csv and Path(lefse_csv).exists():
        lefse_root = Path(lefse_csv).parent.resolve()

    overall_html = ''
    pair_html = ''

    if lefse_root:
        # Overall
        bundle = []
        for c in [lefse_root / "lefse_de_output.csv",
                  lefse_root / "lefse_bar_genus.png",
                  lefse_root / "lefse_dot_genus.png"]:
            if c.exists(): bundle.append(str(c))
        if lefse_csv and Path(lefse_csv).exists() and str(lefse_csv) not in bundle:
            bundle.append(str(Path(lefse_csv)))
        mapped_overall = materialize(bundle, assets / "lefse/overall", mode) if bundle else []
        if mapped_overall:
            imgs = [p for p in mapped_overall if str(p).lower().endswith(IMG_EXT)]
            data = [p for p in mapped_overall if str(p).lower().endswith((".csv",".tsv",".txt"))]
            overall_html = (
                gallery("Overall (Genus)", page_dir_lefse, rels_from(page_dir_lefse, imgs), with_source_links=False) +
                data_card("Overall table", page_dir_lefse, data)
            )

        # Pairwise
        pair_cards = []
        for csvp in sorted(lefse_root.glob("lefse_*vs*.csv")):
            suffix = csvp.stem.replace("lefse_", "")
            bar_png = lefse_root / f"lefse_bar_genus_{suffix}.png"
            dot_png = lefse_root / f"lefse_dot_genus_{suffix}.png"
            bundle = [str(csvp)]
            if bar_png.exists(): bundle.append(str(bar_png))
            if dot_png.exists(): bundle.append(str(dot_png))
            mapped = materialize(bundle, assets / f"lefse/pairwise/{suffix}", mode)
            imgs = [p for p in mapped if str(p).lower().endswith(IMG_EXT)]
            data = [p for p in mapped if str(p).lower().endswith((".csv",".tsv",".txt"))]
            card = gallery(f"Pairwise • {suffix}", page_dir_lefse, rels_from(page_dir_lefse, imgs), with_source_links=False) + \
                   data_card(f"{suffix} table", page_dir_lefse, data)
            pair_cards.append(card)
        if pair_cards:
            pair_html = "<hr/>\n" + "\n".join(pair_cards)

        # 부족 시 스캔 보완
        if not overall_html or not pair_html:
            scan = scan_lefse_images(lefse_root)
            if not overall_html and (scan["overall_bar"] or scan["overall_dot"]):
                mbar = materialize(scan["overall_bar"], assets / "lefse/overall/bar", mode)
                mdot = materialize(scan["overall_dot"], assets / "lefse/overall/dot", mode)
                hb = rels_from(page_dir_lefse, mbar)
                hd = rels_from(page_dir_lefse, mdot)
                parts = ['<section class="section"><h2>Overall (scanned)</h2><div class="card">']
                if hb: parts.append(gallery("Barplot", page_dir_lefse, hb, with_source_links=False))
                if hd: parts.append(gallery("Dotplot", page_dir_lefse, hd, with_source_links=False))
                parts.append('</div></section>')
                overall_html += "\n".join(parts)

            if not pair_html and scan["pairwise"]:
                parts = ['<section class="section"><h2>Pairwise (scanned)</h2>']
                for comp, dd in scan["pairwise"].items():
                    mbar = materialize(dd.get("bar", []), assets / f"lefse/pairwise/{comp}/bar", mode)
                    mdot = materialize(dd.get("dot", []), assets / f"lefse/pairwise/{comp}/dot", mode)
                    parts.append('<div class="card">')
                    parts.append(f'<div class="badge">{comp}</div>')
                    if mbar:
                        parts.append(gallery("Barplot", page_dir_lefse, rels_from(page_dir_lefse, mbar), with_source_links=False))
                    if mdot:
                        parts.append(gallery("Dotplot", page_dir_lefse, rels_from(page_dir_lefse, mdot), with_source_links=False))
                    parts.append('</div>')
                parts.append('</section>')
                pair_html += "\n" + "\n".join(parts)

    if not overall_html:
        overall_html = '<section class="section"><h2>Overall</h2><div class="note">No overall LEfSe outputs found.</div></section>'
    if not pair_html:
        pair_html = '<section class="section"><h2>Pairwise</h2><div class="note">No pairwise LEfSe outputs found.</div></section>'

    write_html(outdir/"lefse/index.html", "LEfSe", root_href="../",
               breadcrumbs_html='Home / LEfSe',
               inner_html=overall_html + "\n" + pair_html)

    # ---------- PICRUSt2 (PDF=errorbar_by_class.pdf만) ----------
    page_dir = outdir / "picrust"
    scan = scan_picrust_dir(Path(pic_plot_dir))
    mapped_figs = materialize(scan["figs"], assets / "picrust", mode) if scan["figs"] else []
    mapped_data = materialize(scan["data"], assets / "picrust", mode) if scan["data"] else []

    hrefs = rels_from(page_dir, mapped_figs)
    inner = gallery("errorbar_by_class", page_dir, hrefs, with_source_links=False)
    inner += data_card("Data sources (KO→pathway & DAA)", page_dir, mapped_data)
    write_html(outdir/"picrust/index.html", "PICRUSt2 Pathway", root_href="../",
               breadcrumbs_html='Home / PICRUSt2 Pathway',
               inner_html=inner)

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ma_outdir", required=True, help="Root of MicrobiomeAnalyst output (contains 01_Alpha, 02_Beta, 03_Composition/props)")
    ap.add_argument("--ma_fig_dir", required=True, help="Flattened images folder (e.g., MA_OUTDIR/99_Figures)")
    ap.add_argument("--lefse_csv", required=False, default=None, help="Path to LEfSe CSV (e.g., MA_OUTDIR/04_LEfse/lefse_de_output.csv)")
    ap.add_argument("--lefse_dir", required=False, default=None, help="Directory containing LEfSe outputs (overall & pairwise)")
    ap.add_argument("--picrust_plot_dir", required=True, help="Folder containing plots_by_class/errorbar_by_class.pdf and TSVs")
    ap.add_argument("--outdir", required=True, help="Output report directory")
    ap.add_argument("--mode", choices=["copy","hardlink"], default="copy")
    args = ap.parse_args()

    build_site(
        ma_outdir=args.ma_outdir,
        ma_fig_dir=args.ma_fig_dir,
        lefse_csv=args.lefse_csv,
        lefse_dir=args.lefse_dir,
        pic_plot_dir=args.picrust_plot_dir,
        outdir=args.outdir,
        mode=args.mode
    )
    print(f"[OK] Built hierarchical report at {Path(args.outdir).resolve()}/index.html")

if __name__ == "__main__":
    main()

