#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:?Usage: $0 <OUTDIR> <TGT>}"
TGT="${2:?Usage: $0 <OUTDIR> <TGT>}"

mkdir -p "$TGT"

echo "# === 1) NSTI 통계 ==="
# gawk로 요약(평균/중앙값/사분위/최소/최대)
gzip -dc "$OUTDIR/marker_predicted_and_nsti.tsv.gz" | \
  awk 'NR>1{v[NR-1]=$NF; n++; s+=$NF; if(min==""||$NF<min)min=$NF; if($NF>max)max=$NF}
       END{
         # 정렬
         asort(v)
         # 중앙값 + 사분위 (단순 분위수 추정)
         if(n%2){med=v[(n+1)/2]} else {med=(v[n/2]+v[n/2+1])/2}
         q1  = v[int((n+3)/4)]
         q3  = v[int((3*n+1)/4)]
         printf("n=%d mean=%.6f median=%.6f q1=%.6f q3=%.6f min=%.6f max=%.6f\n", n, s/n, med, q1, q3, min, max)
       }' | tee "$TGT/nsti_summary.txt"

echo "# === 2) 헤더 확인 ==="
echo -n "KO_unstrat header: "
gzip -dc "$OUTDIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" | head -n 1 || true
echo -n "Pathway_unstrat header: "
gzip -dc "$OUTDIR/pathways_out/path_abun_unstrat.tsv.gz" | head -n 1 || true

echo "# === 3) export (gunzip→TSV, unstrat only) ==="
# KO
gzip -dc "$OUTDIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" > "$TGT/KO_unstrat.tsv"
# EC (있으면)
if [[ -f "$OUTDIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" ]]; then
  gzip -dc "$OUTDIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" > "$TGT/EC_unstrat.tsv"
fi
# Pathways
gzip -dc "$OUTDIR/pathways_out/path_abun_unstrat.tsv.gz" > "$TGT/pathway_unstrat.tsv"

echo "# Done: QC + export (unstrat only)"

