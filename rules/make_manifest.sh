#!/usr/bin/env bash
set -euo pipefail

# 입력(환경변수)
ANALYSIS="${ANALYSIS:-blast_results}"
FASTQ_DIR="${FASTQ_DIR:-/path/to/NCBI_fastqs}"
SUFFIX_R1="${SUFFIX_R1:-_L001_R1_001.fastq.gz}"   # 예: "_1.fastq.gz"
SUFFIX_R2="${SUFFIX_R2:-_L001_R2_001.fastq.gz}"   # 예: "_2.fastq.gz"
MANIFEST="${MANIFEST:-$ANALYSIS/sample_manifest.txt}"

mkdir -p "$(dirname "$MANIFEST")"

echo "[manifest] FASTQ_DIR=$FASTQ_DIR"
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

# R1을 기준으로 짝 찾기 → sample-id는 파일명에서 SUFFIX_R1 제거한 값
found_any=0
while IFS= read -r -d '' r1; do
  r2="${r1%$SUFFIX_R1}$SUFFIX_R2"
  if [[ ! -f "$r2" ]]; then
    echo "[warn] R2 not found for: $r1" >&2
    continue
  fi
  sid="$(basename "$r1")"
  sid="${sid%$SUFFIX_R1}"    # 예: C01_S10_L001_R1_001.fastq.gz → C01_S10_L001  (필요 시 L001 제거 로직 추가 가능)
  # 만약 sample-id를 C01_S10 형태로만 원한다면, 아래 한 줄을 켜세요:
  # sid="${sid%%_L*}"   # 첫 '_L' 이후를 잘라냄 → C01_S10

  printf "%s\t%s\t%s\n" \
    "$sid" \
    "$(readlink -f "$r1")" \
    "$(readlink -f "$r2")" >> "$MANIFEST"
  found_any=1
done < <(find "$FASTQ_DIR" -type f -name "*$SUFFIX_R1" -print0 | sort -z)

if [[ $found_any -eq 0 ]]; then
  echo "ERROR: No R1/R2 pairs found under $FASTQ_DIR (pattern: *$SUFFIX_R1 / *$SUFFIX_R2)" >&2
  exit 1
fi

echo "[manifest] wrote: $MANIFEST"

