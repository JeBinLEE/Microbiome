#!/usr/bin/env bash
set -euo pipefail

ANALYSIS="${ANALYSIS:-vsearch_results}"
REF_DIR="${REF_DIR:-NCBI-RefSeq-16s-202505}"
THREADS_CUTADAPT="${THREADS_CUTADAPT:-24}"
THREADS_BLAST="${THREADS_BLAST:-48}"
MANIFEST="${MANIFEST:-$ANALYSIS/sample_manifest.txt}"

# ---- 경로들 ----
DEMUX="$ANALYSIS/paired-end-demux.qza"
TRIM_QZA="$ANALYSIS/trimmed_sequences.qza"
TRIM_DIR="$ANALYSIS/trimmed_sequences"
STATS_QZA="$ANALYSIS/stats-dada2-2.qza"
TABLE_QZA="$ANALYSIS/table-dada2-2.qza"
REP_QZA="$ANALYSIS/rep-seqs-dada2-2.qza"

# ---- 매니페스트 존재 검사 ----
[[ -s "$MANIFEST" ]] || { echo "ERROR: missing manifest: $MANIFEST"; exit 1; }

# ---- Import ----
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "$DEMUX" \
  --input-format PairedEndFastqManifestPhred33V2

# ---- Cutadapt ----
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$DEMUX" \
  --p-cores "$THREADS_CUTADAPT" \
  --p-front-f "CCTACGGGNGGCWGCAG" \
  --p-front-r "GACTACHVGGGTATCTAATCC" \
  --o-trimmed-sequences "$TRIM_QZA"

# ---- Demux summarize + export ----
mkdir -p "$TRIM_DIR"
qiime demux summarize --i-data "$TRIM_QZA" --o-visualization "$TRIM_DIR/trimmed_sequences.qzv"
qiime tools export --input-path "$TRIM_DIR/trimmed_sequences.qzv" --output-path "$TRIM_DIR"

# ---- DADA2 ----
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$TRIM_QZA" \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 260 --p-trunc-len-r 180 \
  --p-n-threads 0 --p-max-ee-r 5 \
  --o-representative-sequences "$REP_QZA" \
  --o-table "$TABLE_QZA" \
  --o-denoising-stats "$STATS_QZA"

# ---- Export 요약물 ----
STATS_DIR="$ANALYSIS/statsdada2";   mkdir -p "$STATS_DIR"
qiime metadata tabulate --m-input-file "$STATS_QZA" --o-visualization "$STATS_DIR/statsdada2-2.qzv"
qiime tools export --input-path "$STATS_DIR/statsdada2-2.qzv" --output-path "$STATS_DIR"

TABLE_DIR="$ANALYSIS/tabledada2";   mkdir -p "$TABLE_DIR"
qiime feature-table summarize --i-table "$TABLE_QZA" --o-visualization "$TABLE_DIR/tabledada2-2.qzv"
qiime tools export --input-path "$TABLE_DIR/tabledada2-2.qzv" --output-path "$TABLE_DIR"

REP_DIR="$ANALYSIS/rep-seqs-dada2"; mkdir -p "$REP_DIR"
qiime feature-table tabulate-seqs --i-data "$REP_QZA" --o-visualization "$REP_DIR/rep-seqs-dada2-2.qzv"
qiime tools export --input-path "$REP_DIR/rep-seqs-dada2-2.qzv" --output-path "$REP_DIR"

# ---- BLAST 분류 1차 ----
BLASTDB="$REF_DIR/ncbi-refseqs-derep.qza"
TAX_DEREP="$REF_DIR/ncbi-refseqs-taxonomy-derep.qza"
SEARCH1="$ANALYSIS/vsearch.qza"
TAX1="$ANALYSIS/taxonomy_vsearch.qza"
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$REP_QZA" \
  --i-reference-reads "$BLASTDB" \
  --i-reference-taxonomy "$TAX_DEREP" \
  --p-perc-identity 0.85 --p-maxaccepts 1 --p-query-cov 0.85 \
  --p-threads "$THREADS_BLAST" \
  --o-search-results "$SEARCH1" \
  --o-classification "$TAX1"

# ---- 컨타민 제외 & seq 필터 ----
FILT_TABLE="$ANALYSIS/filt_contam_table-dada2_blast.qza"
qiime taxa filter-table \
  --i-table "$TABLE_QZA" \
  --i-taxonomy "$TAX1" \
  --p-exclude mitochondria,chloroplast,Unassigned \
  --o-filtered-table "$FILT_TABLE"

FILT_REP="$ANALYSIS/filtered_rep-seqs-dada2_blast.qza"
qiime feature-table filter-seqs \
  --i-data "$REP_QZA" \
  --i-table "$FILT_TABLE" \
  --o-filtered-data "$FILT_REP"

# ---- BLAST 분류 최종 ----
SEARCH2="$ANALYSIS/vsearch_result.qza"
TAX_FINAL="$ANALYSIS/final_taxonomy_vsearch.qza"
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$FILT_REP" \
  --i-reference-reads "$BLASTDB" \
  --i-reference-taxonomy "$TAX_DEREP" \
  --p-perc-identity 0.85 --p-maxaccepts 1 --p-query-cov 0.85 \
  --p-threads "$THREADS_BLAST" \
  --o-search-results "$SEARCH2" \
  --o-classification "$TAX_FINAL"

# ---- 결과 Export ----
FINAL_DIR="$ANALYSIS/final_taxonomy_vsearch"; mkdir -p "$FINAL_DIR"
qiime metadata tabulate --m-input-file "$TAX_FINAL" --o-visualization "$FINAL_DIR/final_taxonomy_vsearch.qzv"
qiime tools export --input-path "$FINAL_DIR/final_taxonomy_vsearch.qzv" --output-path "$FINAL_DIR"

BR2_DIR="$ANALYSIS/vsearch_results_2"; mkdir -p "$BR2_DIR"
qiime metadata tabulate --m-input-file "$SEARCH2" --o-visualization "$BR2_DIR/vsearch_results_2.qzv"
qiime tools export --input-path "$BR2_DIR/vsearch_results_2.qzv" --output-path "$BR2_DIR"

ASV_DIR="$ANALYSIS/ASV_quantified";  mkdir -p "$ASV_DIR"
qiime tools export --input-path "$FILT_TABLE" --output-path "$ASV_DIR"
biom convert -i "$ASV_DIR/feature-table.biom" -o "$ASV_DIR/feature-table.txt" --to-tsv

TB_DIR="$ANALYSIS/taxonomy_barplot"; mkdir -p "$TB_DIR"
qiime taxa barplot --i-table "$FILT_TABLE" --i-taxonomy "$TAX_FINAL" --o-visualization "$TB_DIR/taxonomy_barplot.qzv"
qiime tools export --input-path "$TB_DIR/taxonomy_barplot.qzv" --output-path "$TB_DIR"

echo "[analysis] done"

