# Updated microbiome Snakemake files (sweeps + filtering + norm)

Included:
- `Snakefile3` (updated): runs MicrobiomeAnalyst-style pipeline **per normalization sweep** and publishes a **primary** run to the original output location to keep downstream rules (figures/report) working.
- `config.yaml` (updated): adds `filter`, `alpha`, `beta`, `sweeps`, `lefse` sections and points `script_ma_pipeline_r` to the new R script under `scripts/`.
- `scripts/microbiomeanalyst_pipeline.R` (new): filtering + multiple normalization options, alpha/beta analyses, and composition outputs.
- `scripts/run_lefse.py` (copied): optional LEfSe sweeps (python runner) reading filtered tables produced by each normalization run.

Key paths:
- per-norm outputs: `${MA_OUTDIR}/sweeps/<norm_id>/`
- primary published outputs: `${MA_OUTDIR}/00_QC_Normalization`, `01_Alpha`, `02_Beta`, `03_Composition`

How to run:
- `snakemake -s Snakefile3 -j <N> --configfile config.yaml`

Configure sweeps:
- edit `sweeps.normalizations` in `config.yaml`
- choose `sweeps.primary_norm` to publish to `${MA_OUTDIR}/`

Configure filtering:
- edit `filter.*` in `config.yaml`

Optional LEfSe sweeps:
- set `lefse.enabled: true` and edit `lefse.sweeps` list in `config.yaml`
- outputs go to `${MA_OUTDIR}/sweeps/<norm_id>/04_LEfse_py/<lefse_id>/`

Changes in v5:
- Legacy `MicrobiomeAnalystR`-based LEfSe (`scripts/run_ma_lefse.R`) is now **optional** and will only run if `lefse.legacy_maR: true` in `config.yaml`.
- `scripts/microbiomeanalyst_pipeline.R` no longer depends on `broom` (avoids missing-package crashes).
- Filtering summary now includes `00_QC_Normalization/filter_details.tsv` and `qc_summary_long.tsv` so you can directly report how many features remain after each filtering step.
