#!/usr/bin/env bash
set -euo pipefail

REF_DIR="${REF_DIR:-NCBI-RefSeq-16s-202505}"
mkdir -p "$REF_DIR"

BLASTDB="$REF_DIR/ncbi-refseqs-blastdb.qza"
DEREP_TAX="$REF_DIR/ncbi-refseqs-taxonomy-derep.qza"

# 최종 산출물 있으면 즉시 스킵
if [[ -f "$BLASTDB" && -f "$DEREP_TAX" ]]; then
  echo "[refdb] Final outputs exist → skip rebuild"
  exit 0
fi

SEQS_UNF="$REF_DIR/ncbi-refseqs-unfiltered.qza"
TAX_UNF="$REF_DIR/ncbi-refseqs-taxonomy-unfiltered.qza"
KEPT="$REF_DIR/ncbi-refseqs.qza"
DROP="$REF_DIR/ncbi-refseqs-too-short.qza"
TAX_CLEAN="$REF_DIR/ncbi-refseqs-taxonomy-clean.qza"
DEREP_SEQS="$REF_DIR/ncbi-refseqs-derep.qza"

# 한 개의 인자로 전달되도록 큰따옴표로 감쌈
QUERY='(33175[BioProject]) OR (33317[BioProject])'

echo "[refdb] get-ncbi-data..."
qiime rescript get-ncbi-data \
  --p-query "$QUERY" \
  --o-sequences "$SEQS_UNF" \
  --o-taxonomy "$TAX_UNF"

echo "[refdb] filter-seqs-length-by-taxon..."
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences "$SEQS_UNF" \
  --i-taxonomy "$TAX_UNF" \
  --p-labels Bacteria Archaea \
  --p-min-lens 1200 1200 \
  --o-filtered-seqs "$KEPT" \
  --o-discarded-seqs "$DROP"

echo "[refdb] taxonomy cleanup..."
qiime rescript filter-taxa \
  --i-taxonomy "$TAX_UNF" \
  --p-exclude 'uncultured,unclassified,metagenome,environmental,sample,unknown,incertae sedis' \
  --m-ids-to-keep-file "$KEPT" \
  --o-filtered-taxonomy "$TAX_CLEAN"

echo "[refdb] dereplicate..."
qiime rescript dereplicate \
  --i-sequences "$KEPT" \
  --i-taxa "$TAX_CLEAN" \
  --p-mode 'uniq' \
  --o-dereplicated-sequences "$DEREP_SEQS" \
  --o-dereplicated-taxa "$DEREP_TAX"

echo "[refdb] makeblastdb..."
qiime feature-classifier makeblastdb \
  --i-sequences "$DEREP_SEQS" \
  --o-database "$BLASTDB"

echo "[refdb] done → $BLASTDB  /  $DEREP_TAX"

