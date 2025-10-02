# MicrobiomeAnalystR + QIIME2 Snakemake Pipeline

End-to-end workflow for 16S microbiome analysis:

1) QIIME2 preprocessing  
2) MicrobiomeAnalystR “full” analysis  
3) Pairwise LEfSe (ALL + every group pair)  
4) Figure collection & flattening

This repository tracks only the workflow code and essential artifacts (e.g., YAML env, `rules/`, `scripts/`, `Snakefile_refine`). Result data are intentionally not versioned.

---

## Repository Contents

- `Snakefile_refine` – main Snakemake workflow  
- `rules/` – shell helpers (`2_NCBI_build_db.sh`, `3_analyze.sh`, `make_manifest.sh`, `make_MicroAnalysis_input.sh`)  
- `scripts/` – R scripts (`microbiomeanalyst_pipeline.R`, `run_ma_lefse.R`)  
- `environment.yml` – (optional) conda/mamba environment spec (rename to your actual file if different)

---

## Dependencies

- **Snakemake**
- **QIIME2**
- **R (Rscript)** with: `MicrobiomeAnalystR`, `phyloseq`, `dplyr`, `ggplot2`, `readr`

---

## 1) Install via Conda/Mamba (YAML)

If you maintain an environment file (e.g., `environment.yml`):

```bash
# faster with mamba
mamba env create -f qiime2-amplicon-ubuntu-latest-conda.yml
# or with conda
conda env create -f qiime2-amplicon-ubuntu-latest-conda.yml
```
## 1-2) Installing MicrobiomeAnalystR (R)

MicrobiomeAnalystR depends on several R/Bioconductor packages; R and Bioconductor version compatibility matters. The project and the following issue contain useful installation guidance:

Project: https://github.com/xia-lab/MicrobiomeAnalystR

Issue discussion (worth reading): https://github.com/xia-lab/MicrobiomeAnalystR/issues/35

## Environment Variables
You can override these before running Snakemake:

| Variable           | Default                          | Description                                                      |
| ------------------ | -------------------------------- | ---------------------------------------------------------------- |
| `REF_DIR`          | `NCBI-RefSeq-16s-202505`         | Path to the 16S reference database                               |
| `ANALYSIS`         | `GutBiomeTech/blast_results`     | Root directory for QIIME2 analysis outputs                       |
| `FASTQ_DIR`        | `/home/ljb/.../M17_fastq`        | Location of FASTQ files                                          |
| `SUFFIX_R1`        | `_L001_R1_001.fastq.gz`          | R1 filename suffix                                               |
| `SUFFIX_R2`        | `_L001_R2_001.fastq.gz`          | R2 filename suffix                                               |
| `THREADS_CUTADAPT` | `24`                             | Threads for cutadapt                                             |
| `THREADS_BLAST`    | `48`                             | Threads for BLAST                                                |
| `MA_DIR`           | `GutBiomeTech/MicrobiomeAnalyst` | Base path for MicrobiomeAnalyst files                            |
| `MA_OUTDIR`        | `${MA_DIR}/results_R_test_MA2`   | Output dir for the full MicrobiomeAnalyst pipeline               |
| `MA_GROUP`         | `Group1`                         | Metadata column used as experimental group                       |
| `LEFSE_OUTDIR`     | `/home/.../04_LEfse`             | Output dir for pairwise LEfSe results                            |
| `RS_BIN`           | `/usr/bin/Rscript`               | Path to Rscript binary                                           |

## 2) How to Run
```bash
snakemake -j 8 --snakefile Snakemake_refine
```


## 3) What Each Rule Does (Snakefile)

### rule all
- Collects the final targets:

- QIIME2 pipeline sentinel: .done/_3_analyze_sh.done

- MicrobiomeAnalyst upload files: upload/feature_table_for_MA.txt, metadata_for_MA.txt, taxonomy_for_MA.txt

- Full MicrobiomeAnalyst pipeline sentinel: ${MA_OUTDIR}/.done_pipeline

- Pairwise LEfSe sentinel: ${LEFSE_OUTDIR}/.done_pairwise

- Figure collection sentinel: ${MA_OUTDIR}/.done_figures

### 1) build_db
- Builds the NCBI RefSeq 16S database artifacts. Skipped if outputs already exist.
- Outputs:

- ${REF_DIR}/ncbi-refseqs-blastdb.qza

- ${REF_DIR}/ncbi-refseqs-taxonomy-derep.qza

### 1.5) make_manifest
- Generates sample_manifest.txt from the FASTQ directory, using SUFFIX_R1/SUFFIX_R2.
- Input depends on build_db (ordering only).

### 2) analyze
- Runs cutadapt/BLAST etc. via rules/3_analyze.sh using the manifest and reference DB.
- Output sentinel: .done/_3_analyze_sh.done.

### 3) make_ma_upload_from_manifest_sh
- Creates three MicrobiomeAnalyst input tables:

1. feature_table_for_MA.txt

2. taxonomy_for_MA.txt

3. metadata_for_MA.txt
- from QIIME2 outputs using rules/make_MicroAnalysis_input.sh.

### 4) run_microbiomeanalyst_pipeline
- Runs full MicrobiomeAnalystR analyses on all samples via scripts/microbiomeanalyst_pipeline.R.
- Parameters: --group (e.g., Group1), --rarefy_q (min or quantile).
- Output sentinel: ${MA_OUTDIR}/.done_pipeline.

### 5) run_microbiomeanalyst_pairwise
- Runs ALL + pairwise LEfSe to a dedicated folder via scripts/run_ma_lefse.R.
- Parameters: --group, --rarefy_q (getopt-based).
- Output sentinel: ${LEFSE_OUTDIR}/.done_pairwise.
- Typical artifacts: lefse_ALL_Genus.csv, lefse_bar_genus_C_vs_D.png, etc.

### 6) collect_ma_figures
- Flattens all *.png/*.pdf under ${MA_OUTDIR} into ${MA_OUTDIR}/99_Figures/, and creates _index.tsv.
- Output sentinel: ${MA_OUTDIR}/.done_figures.

