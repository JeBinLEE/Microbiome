# 🧬 Microbiome Analysis Pipeline

A Snakemake-based microbiome analysis pipeline with **normalization sweeps** for 16S rRNA amplicon sequencing data.

---

## 📋 Overview

This pipeline performs comprehensive microbiome analysis including:
- QIIME2-based sequence processing
- Multiple normalization strategies for different analysis types
- Alpha/Beta diversity analysis
- Composition analysis with group-level visualizations
- LEfSe biomarker discovery
- PICRUSt2 functional prediction

---

## 🗂️ Repository Structure

```
├── Snakefile3                    # Main Snakemake pipeline
├── config.yaml                   # Configuration file
├── envs/                         # Conda environment specs
│   └── qiime2-amplicon-*.yml
├── rules/                        # Shell helper scripts
│   ├── 2_NCBI_build_db.sh
│   ├── 3_analyze.sh
│   ├── make_manifest.sh
│   └── make_MicroAnalysis_input.sh
└── scripts/                      # Analysis scripts
    ├── microbiomeanalyst_pipeline.R   # Core analysis (R)
    ├── run_lefse.py                   # LEfSe analysis (Python)
    ├── run_ma_lefse.R                 # LEfSe via MicrobiomeAnalystR
    ├── plot_kegg_by_class.R           # KEGG pathway plots
    └── build_html_report.py           # HTML report generator
```

---

## 🔧 Normalization Sweeps

Different analyses require different normalization methods. This pipeline supports multiple normalizations in a single run:

| Sweep ID | Method | Use Case |
|----------|--------|----------|
| `01_AlphaBeta_Rarefy40k` | Rarefaction (40,000 reads) | Alpha/Beta diversity |
| `02_Composition_TSS` | TSS (Total Sum Scaling) | Stacked bar plots, LEfSe |
| `03_Heatmap_CLR` | CLR (Centered Log-Ratio) | Heatmap, Correlation |
| `04_PICRUSt_Raw` | Raw counts | PICRUSt2 functional prediction |

### Output Structure

```
results/sweeps/
├── 01_AlphaBeta_Rarefy40k/
│   ├── 00_QC_Normalization/
│   ├── 01_Alpha/           # Alpha diversity results
│   ├── 02_Beta/            # Beta diversity results
│   └── 03_Composition/
├── 02_Composition_TSS/
│   └── 03_Composition/     # Group-averaged bar plots
│       ├── stacked_bar_by_sample.pdf
│       ├── stacked_bar_by_group.pdf    # ✨ Group average
│       ├── area_plot_by_group.pdf      # ✨ Area plot
│       └── composition_by_group.tsv
├── 03_Heatmap_CLR/
│   └── 00_QC_Normalization/
│       └── feature_table_normalized.tsv  # CLR-transformed data
└── 04_PICRUSt_Raw/
```

---

## 🚀 Quick Start

### 1. Install Dependencies

```bash
# Create conda environment
mamba env create -f envs/qiime2-amplicon-ubuntu-latest-conda.yml
conda activate qiime2-re

# Install R packages (if needed)
Rscript -e "install.packages(c('optparse', 'phyloseq', 'vegan', 'ggplot2', 'data.table'))"
```

### 2. Configure

Edit `config.yaml`:

```yaml
# Input/Output paths
ma_dir: "path/to/MicrobiomeAnalyst"
ma_outdir: "path/to/results"
ma_group: "Group1"  # Metadata column for grouping

# Normalization sweeps
sweeps:
  primary_norm: "01_AlphaBeta_Rarefy40k"
  normalizations:
    - id: "01_AlphaBeta_Rarefy40k"
      scale: "rarefy"
      rarefy_depth: 40000
      # ... more options
```

### 3. Run Pipeline

```bash
snakemake -j 40 --snakefile Snakefile3 --configfile config.yaml
```

---

## 📊 Analysis Outputs

### Alpha Diversity
- Observed richness, Chao1, Shannon, Simpson indices
- Boxplots with statistical tests (Kruskal-Wallis, ANOVA)
- Location: `sweeps/01_AlphaBeta_Rarefy40k/01_Alpha/`

### Beta Diversity  
- Distance matrices: Bray-Curtis, Jaccard, Euclidean
- Ordination: PCoA, NMDS
- Statistical tests: PERMANOVA, ANOSIM
- Location: `sweeps/01_AlphaBeta_Rarefy40k/02_Beta/`

### Composition
- Stacked bar plots (sample-level and group-averaged)
- Area plots for composition trends
- Top 20 taxa visualization
- Location: `sweeps/02_Composition_TSS/03_Composition/`

### LEfSe
- Pairwise group comparisons
- LDA scores for biomarker identification
- Location: `results/04_LEfse/`

### PICRUSt2
- KO (KEGG Ortholog) predictions
- Pathway abundance tables
- KEGG class visualizations
- Location: `picrust2/`

---

## ⚙️ Configuration Options

### Filtering Parameters

```yaml
filter:
  min_total_count: 4      # Min total reads per feature
  prevalence: 0.1         # Min sample prevalence (10%)
  iqr_remove: 0.1         # Remove low-variance features
```

### Alpha Diversity

```yaml
alpha:
  ranks: "ASV,Genus,Family,Species"
  measures: "Observed,Chao1,Shannon,Simpson"
  stats: "kruskal,anova,wilcoxon"
```

### Beta Diversity

```yaml
beta:
  ranks: "Genus,Family"
  distances: "bray,jaccard,euclidean"
  ordinations: "PCoA,NMDS"
  stats: "permanova,anosim"
  permutations: 999
```

---

## 📦 Dependencies

### Required Software
- Snakemake ≥7.0
- QIIME2 (2024.5+)
- R ≥4.0
- Python ≥3.8

### R Packages
- `phyloseq`, `vegan`, `ggplot2`, `dplyr`, `tidyr`
- `data.table`, `optparse`, `pheatmap`
- `MicrobiomeAnalystR` (optional)

### Python Packages
- `pandas`, `numpy`, `biom-format`
- `ggpicrust2` (for KEGG visualization)

---

## 📝 Notes

### TSS/CLR and Alpha Diversity
Alpha diversity metrics (Observed, Chao1) require **integer count data**. When using TSS or CLR normalization, alpha diversity calculation is automatically skipped with a note file generated.

### Customizing Sweeps
Add new normalization sweeps in `config.yaml`:

```yaml
sweeps:
  normalizations:
    - id: "custom_rarefy_10k"
      scale: "rarefy"
      rarefy_depth: 10000
      transform: "none"
```

---

## 📄 License

MIT License

## 👤 Author

JeBin Lee - [GitHub](https://github.com/JeBinLEE)
