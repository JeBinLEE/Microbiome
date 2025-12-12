## Repository Contents

- `Snakefile` – main Snakemake workflow
- `config.yaml` – configuration
- `envs/` – conda environment specs
- `rules/` – shell helpers  
  - `2_NCBI_build_db.sh`  
  - `3_analyze.sh`  
  - `make_manifest.sh`  
  - `make_MicroAnalysis_input.sh`
- `scripts/` – R scripts  
  - `microbiomeanalyst_pipeline.R`  
  - `run_ma_lefse.R`
- `pipeline.pdf` – workflow diagram

> This repository tracks only the workflow code and essential artifacts.  
> Result data are intentionally not versioned.

---

## Dependencies

- Snakemake
- QIIME2
- R (Rscript)

R packages used in this pipeline:
`readr`, `dplyr`, `tibble`, `stringr`, `ggplot2`, `ggprism`, `ggpicrust2`, `rlang`,  
`optparse`, `data.table`, `pheatmap`, `reshape2`, `vegan`, `phyloseq`, `withr`, `tidyr`,  
`MicrobiomeAnalystR`.

### Note: Installing Tax4Fun on Linux

`Tax4Fun` is no longer available on CRAN, so it needs to be installed from the original website as a source tarball.


1. Download the tarball:

```bash
wget http://tax4fun.gobics.de/Tax4Fun/Tax4Fun_0.3.1.tar.gz
install.packages("Tax4Fun_0.3.1.tar.gz",
                 repos = NULL,
                 type  = "source")
```

### Note: Installing MicrobiomeAnalystR
https://github.com/xia-lab/MicrobiomeAnalystR/issues/35
---

## 1) Install via Conda/Mamba

Use the environment file under `envs/`:

```bash
# faster with mamba
mamba env create -f envs/environment.yml

# or with conda
conda env create -f envs/environment.yml
```


```bash
conda activate <ENV_NAME>

```
## 1-2) Installing MicrobiomeAnalystR (R)

MicrobiomeAnalystR depends on several R/Bioconductor packages; R and Bioconductor version compatibility matters. The project and the following issue contain useful installation guidance:

Project: https://github.com/xia-lab/MicrobiomeAnalystR

Issue discussion (worth reading): https://github.com/xia-lab/MicrobiomeAnalystR/issues/35

## Environment Variables
You can override these before running Snakemake:

| Variable           | Default                                                                                | Description                                        |
| ------------------ | -------------------------------------------------------------------------------------- | -------------------------------------------------- |
| `REF_DIR`          | `NCBI-RefSeq-16s-202505`                                                               | Path/name of the 16S reference database            |
| `ANALYSIS`         | `GutBiomeTech/vsearch_results`                                                         | Root directory for QIIME2 analysis outputs         |
| `FASTQ_DIR`        | `/home/ljb/qiime2_analysis/GutBiomeTech/M17_fastq`                                     | Location of FASTQ files                            |
| `SUFFIX_R1`        | `_L001_R1_001.fastq.gz`                                                                | R1 filename suffix                                 |
| `SUFFIX_R2`        | `_L001_R2_001.fastq.gz`                                                                | R2 filename suffix                                 |
| `THREADS_CUTADAPT` | `24`                                                                                   | Threads for cutadapt                               |
| `THREADS_BLAST`    | `48`                                                                                   | Threads for BLAST                                  |
| `MA_DIR`           | `GutBiomeTech/MicrobiomeAnalyst`                                                       | Base path for MicrobiomeAnalyst files              |
| `MA_OUTDIR`        | `${MA_DIR}/results_R_test_MA2`                                                         | Output dir for the full MicrobiomeAnalyst pipeline |
| `MA_GROUP`         | `Group1`                                                                               | Metadata column used as experimental group         |
| `RAREFY_Q`         | `0.00`                                                                                 | 0.00 = min depth, 0~1 = quantile                   |
| `LEFSE_OUTDIR`     | `/home/ljb/qiime2_analysis/GutBiomeTech/MicrobiomeAnalyst/results_R_test_MA2/04_LEfse` | Output dir for pairwise LEfSe                      |
| `RS_BIN`           | `/usr/bin/Rscript`                                                                     | Path to Rscript binary                             |


## 2) How to Run
```bash
snakemake -j 8 --snakefile Snakefile --config config.yaml
```

## 3) What Each Rule Does (Snakefile)

### rule all
- Collects the final targets:
  - QIIME2 pipeline sentinel:
    - `.done/_3_analyze_sh.done`
  - MicrobiomeAnalyst upload files:
    - `${MA_DIR}/upload/feature_table_for_MA.txt`
    - `${MA_DIR}/upload/metadata_for_MA.txt`
    - `${MA_DIR}/upload/taxonomy_for_MA.txt`
  - Full MicrobiomeAnalyst pipeline sentinel:
    - `${MA_OUTDIR}/.done_pipeline`
  - Pairwise LEfSe sentinel:
    - `${LEFSE_OUTDIR}/.done_pairwise`
  - Figure collection sentinel:
    - `${MA_OUTDIR}/.done_figures`
  - PICRUSt2 core outputs:
    - `${ANALYSIS}/picrust2/.done_picrust2`
    - `${ANALYSIS}/picrust2/picrust2_data/KO_unstrat.tsv`
    - `${ANALYSIS}/picrust2/picrust2_data/pathway_unstrat.tsv`
    - `${ANALYSIS}/picrust2/plots_by_class/errorbar_by_class.pdf`
  - Report:
    - `report/index.html`

---

### 1) build_db
- Builds the NCBI RefSeq 16S database artifacts.
- Skipped if outputs already exist.
- Outputs:
  - `${REF_DIR}/ncbi-refseqs-blastdb.qza`
  - `${REF_DIR}/ncbi-refseqs-taxonomy-derep.qza`
- Runs:
  - `rules/2_NCBI_build_db.sh`

---

### 1.5) make_manifest
- Generates `sample_manifest.txt` from the FASTQ directory.
- Uses `SUFFIX_R1/SUFFIX_R2` to detect paired reads.
- Input depends on `build_db` (ordering only).
- Output:
  - `${ANALYSIS}/sample_manifest.txt`
- Runs:
  - `rules/make_manifest.sh`

---

### 2) run_qiime2
- Runs cutadapt/BLAST and the core QIIME2-based analysis using the manifest and reference DB.
- Key outputs:
  - `.done/_3_analyze_sh.done`
  - `${ANALYSIS}/ASV_quantified/feature-table.txt`
  - `${ANALYSIS}/rep-seqs-dada2-2.qza`
  - `${ANALYSIS}/table-dada2-2.qza`
- Runs:
  - `rules/3_analyze.sh`
- Conda:
  - `envs/qiime2-amplicon-ubuntu-latest-conda.yml`

---

### 3) make_MicrobiomeAnalystR_input
- Creates three MicrobiomeAnalyst upload-ready input files from QIIME2 outputs:
  - `feature_table_for_MA.txt`
  - `taxonomy_for_MA.txt`
  - `metadata_for_MA.txt`
- Outputs:
  - `${MA_DIR}/upload/feature_table_for_MA.txt`
  - `${MA_DIR}/upload/metadata_for_MA.txt`
  - `${MA_DIR}/upload/taxonomy_for_MA.txt`
- Runs:
  - `rules/make_MicroAnalysis_input.sh`

---

### 4) run_MicrobiomeAnalyst_pipeline
- Runs the full MicrobiomeAnalystR pipeline on all samples.
- Script:
  - `scripts/microbiomeanalyst_pipeline.R`
- Parameters:
  - `--group` (e.g., `Group1`)
  - `--rarefy_q` (`0.00` = min depth or `0~1` = quantile)
  - `--p_cut 0.1`
  - `--count_cut 4`
- Output sentinel:
  - `${MA_OUTDIR}/.done_pipeline`

---

### 5) run_MicrobiomeAnalyst_LefSe
- Runs **ALL + pairwise LEfSe** to a dedicated folder.
- Script:
  - `scripts/run_ma_lefse.R`
- Parameters:
  - `--group`
  - `--rarefy_q`
  - `--p_cut 0.1`
  - `--count_cut 4`
- Output sentinel:
  - `${LEFSE_OUTDIR}/.done_pairwise`

---

### 6) collect_ma_figures
- Flattens all `*.png/*.pdf` under `${MA_OUTDIR}` into a single folder and generates an index.
- Outputs:
  - `${MA_OUTDIR}/99_Figures/`
  - `${MA_OUTDIR}/99_Figures/_index.tsv`
- Output sentinel:
  - `${MA_OUTDIR}/.done_figures`

---

## 4) PICRUSt2 Rules

### 7) export_repseqs_fna
- Exports FASTA from QIIME2 rep-seqs.
- Input:
  - `${ANALYSIS}/rep-seqs-dada2-2.qza`
- Output:
  - `${ANALYSIS}/ASV_quantified/rep_seqs.fna`

---

### 8) export_table_biom
- Exports BIOM from QIIME2 feature table.
- Input:
  - `${ANALYSIS}/table-dada2-2.qza`
- Output:
  - `${ANALYSIS}/ASV_quantified/table_seqs.biom`

---

### 9) run_picrust2
- Runs PICRUSt2 prediction pipeline.
- Inputs:
  - `${ANALYSIS}/ASV_quantified/rep_seqs.fna`
  - `${ANALYSIS}/ASV_quantified/table_seqs.biom`
- Outputs:
  - `${ANALYSIS}/picrust2/`
  - `${ANALYSIS}/picrust2/.done_picrust2`
- Conda:
  - `envs/picrust2.yml`

---

### 10) qc_and_export_picrust2
- Performs QC and exports unstratified outputs.
- Outputs:
  - `${ANALYSIS}/picrust2/picrust2_data/KO_unstrat.tsv`
  - `${ANALYSIS}/picrust2/picrust2_data/pathway_unstrat.tsv`
  - `${ANALYSIS}/picrust2/picrust2_data/nsti_summary.txt`
- Script:
  - `scripts/picrust2_qc_and_export.sh`

---

### 11) plot_kegg_by_class
- Generates class-level KEGG errorbar plots using KO results and metadata.
- Inputs:
  - `${ANALYSIS}/picrust2/picrust2_data/KO_unstrat.tsv`
  - `${MA_DIR}/upload/metadata_for_MA.txt`
- Output:
  - `${ANALYSIS}/picrust2/plots_by_class/errorbar_by_class.pdf`
- Script:
  - `scripts/plot_kegg_by_class.R`

---

## 5) Report

### 12) build_html_report
- Builds a consolidated HTML report integrating:
  - MA figure collection
  - Pairwise LEfSe results
  - PICRUSt2 plots and summary tables
- Output:
  - `report/index.html`
- Note:
  - This rule calls a Python script via `SCRIPT_BUILD_REPORT`.
  - Ensure `SCRIPT_BUILD_REPORT` is defined in the Snakefile or exported before running.
