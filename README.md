# Microbiome
microbiome analysis Raw to LEfse

## 의존성
- Snakemake
- QIIME2 (+ cutadapt 등)
- R (Rscript), `MicrobiomeAnalystR`, `phyloseq`, `dplyr`, `ggplot2`, `readr`
- Graphviz(`dot`) – DAG 그릴 때만

## 환경변수
| 변수 | 기본값 | 설명 |
|---|---|---|
| `REF_DIR` | `NCBI-RefSeq-16s-202505` | 16S DB 경로 |
| `ANALYSIS` | `GutBiomeTech/blast_results` | QIIME2 분석 산출 루트 |
| `FASTQ_DIR` | `/home/ljb/.../M17_fastq` | FASTQ 폴더 |
| `SUFFIX_R1` | `_L001_R1_001.fastq.gz` | R1 접미사 |
| `SUFFIX_R2` | `_L001_R2_001.fastq.gz` | R2 접미사 |
| `THREADS_CUTADAPT` | `24` | cutadapt 스레드 |
| `THREADS_BLAST` | `48` | BLAST 스레드 |
| `MA_DIR` | `GutBiomeTech/MicrobiomeAnalyst` | MicrobiomeAnalyst 베이스 |
| `MA_OUTDIR` | `${MA_DIR}/results_R_test_MA2` | MA 전체 파이프라인 출력 |
| `MA_GROUP` | `Group1` | 메타데이터 그룹 열명 |
| `RAREFY_Q` | `0.00` | 0이면 min, 0~1이면 quantile rarefaction |
| `LEFSE_OUTDIR` | `/home/.../04_LEfse` | pairwise LEfSe 출력 디렉토리 |
| `RS_BIN` | `/usr/bin/Rscript` | Rscript 실행파일 |

## 실행 방법

### 0) 환경변수(필요 시)
```bash
export MA_OUTDIR="/path/to/results"
export LEFSE_OUTDIR="/path/to/results/04_LEfse"
export MA_GROUP="Group1"
export RAREFY_Q="0.00"
