#!/bin/sh
# POSIX /bin/sh 호환: manifest 1열만 사용해서 metadata 만들기(+ feature/taxonomy 변환)
set -eu

ANALYSIS="${ANALYSIS:-../GutBiomeTech/vsearch_results}"
MA_DIR="${MA_DIR:-../GutBiomeTech/MicrobiomeAnalyst}"
UP="$MA_DIR/upload"
mkdir -p "$UP"

FT="$ANALYSIS/ASV_quantified/feature-table.txt"
MAN="$ANALYSIS/sample_manifest.txt"
TAX1="$ANALYSIS/final_taxonomy_vsearch/taxonomy.tsv"
TAX2="$ANALYSIS/final_taxonomy_vsearch/metadata.tsv"

[ -s "$FT" ]  || { echo "ERROR: feature-table 없음: $FT"  >&2; exit 1; }
[ -s "$MAN" ] || { echo "ERROR: sample_manifest 없음: $MAN" >&2; exit 1; }
if   [ -s "$TAX1" ]; then TAX_SRC="$TAX1";
elif [ -s "$TAX2" ]; then TAX_SRC="$TAX2";
else echo "ERROR: taxonomy 소스 없음: $TAX1 또는 $TAX2" >&2; exit 1; fi

# 1) feature table → #NAME 규격
#awk 'BEGIN{FS=OFS="\t"} NR==1{$1="#NAME"} {sub("\r$",""); print}' \
#  "$FT" > "$UP/feature_table_for_MA.txt"

awk -F '\t' '
  BEGIN{OFS="\t"; hdr=0}
  {
    sub(/\r$/,"")                       
    if (!hdr && $0 ~ /^[[:space:]]*#/ && $1!="#OTU ID") next
    if (!hdr && ($1=="Feature ID" || $1=="FeatureID" || $1=="#OTU ID" || $1=="OTU ID")) {
      $1="#NAME"; hdr=1; print; next
    }
    if (hdr) { print }
  }
' "$FT" > "$UP/feature_table_for_MA.txt"

# 2) metadata → #NAME, Sample, Group1 (manifest 1열)
# GROUP_CHARS: Number of characters from sample ID to use for grouping (default: 1)
GROUP_CHARS="${GROUP_CHARS:-1}"

{
  echo "#NAME	Sample	Group1"
  awk -v grp_chars="$GROUP_CHARS" '
    BEGIN{ OFS="\t" }
    NR==1 { next }                             # manifest header
    {
      line=$0; sub(/\r$/,"",line)
      n=split(line, a, /[ \t]+/); id=a[1]      
      if (id=="" || tolower(id)=="sample-id") next

      if (length(id) >= grp_chars) {
        grp = substr(id, 1, grp_chars)
      } else {
        grp = id
      }
      if (grp=="") grp="G1"

      print id, id, grp
    }
  ' "$MAN"
} > "$UP/metadata_for_MA.txt"

# (참고용) feature-table의 샘플명 추출: 헤더 라인 인식해서 2열~끝
FT_SAMPLES="$(mktemp)"; trap 'rm -f "$FT_SAMPLES"' EXIT
awk -F '\t' '
  {
    gsub(/\r$/,"")
    if ($1=="Feature ID" || $1=="FeatureID" || $1=="#OTU ID" || $1=="OTU ID") {
      for (i=2; i<=NF; i++) if ($i!="") print $i
      exit
    }
  }
' "$FT" > "$FT_SAMPLES"

## 3) taxonomy → #TAXONOMY + 7 ranks
#awk -F'\t' '
#BEGIN{OFS="\t"}
#/^[[:space:]]*#/ {next}
#NR==1{
#  fi=ti=0
#  for(i=1;i<=NF;i++){ if($i=="Feature ID"||$i=="FeatureID") fi=i }
#  for(i=1;i<=NF;i++){ if($i=="Taxon"||$i=="taxonomy"||$i=="Classification"||$i=="Consensus") ti=i }
#  if(!fi||!ti){print "ERROR: need FeatureID & Taxon columns" > "/dev/stderr"; exit 1}
#  print "#TAXONOMY","Kingdom","Phylum","Class","Order","Family","Genus","Species"
#  next
#}
#{
#  fid=$fi; tax=$ti
#  sub(/^[ \t]+/,"",tax); sub(/[ \t]+$/,"",tax)
#  n=split(tax,a,/;[ ]*/)
#  for(i=1;i<=7;i++) out[i]=""
#  for(i=1;i<=n && i<=7;i++) out[i]=a[i]
#  print fid,out[1],out[2],out[3],out[4],out[5],out[6],out[7]
#}
#' "$TAX_SRC" > "$UP/taxonomy_for_MA.txt"

# 3) taxonomy → #TAXONOMY + (Kingdom 제외) 6 ranks
awk -F'\t' '
BEGIN{OFS="\t"}
/^[[:space:]]*#/ {next}
NR==1{
  fi=ti=0
  for(i=1;i<=NF;i++){ if($i=="Feature ID"||$i=="FeatureID") fi=i }
  for(i=1;i<=NF;i++){ if($i=="Taxon"||$i=="taxonomy"||$i=="Classification"||$i=="Consensus") ti=i }
  if(!fi||!ti){print "ERROR: need FeatureID & Taxon columns" > "/dev/stderr"; exit 1}
  print "#TAXONOMY","Kingdom","Phylum","Class","Order","Family","Genus","Species"
  next
}
{
  fid=$fi; tax=$ti
  sub(/^[ \t]+/,"",tax); sub(/[ \t]+$/,"",tax)
  n=split(tax,a,/;[ ]*/)

  for(i=1;i<=7;i++) out[i]=""
  for(i=1;i<=7;i++){
    if (i<=n) out[i]=a[i]
  }
  print fid,out[1],out[2],out[3],out[4],out[5],out[6],out[7]
}
' "$TAX_SRC" > "$UP/taxonomy_for_MA.txt"


echo "Done. Created files:"
ls -lh "$UP/feature_table_for_MA.txt" "$UP/metadata_for_MA.txt" "$UP/taxonomy_for_MA.txt"

# (선택) 간단 점검
echo "FT sample count: $(wc -l < "$FT_SAMPLES" | tr -d " ")"
echo "MAN sample count: $(awk "NR>1{print 1}" "$MAN" | wc -l | tr -d " ")"

