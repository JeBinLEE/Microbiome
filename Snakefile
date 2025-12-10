from snakemake.shell import shell
shell.executable("/bin/bash")
import os

# =========================================================
# 환경변수(필요 시 export로 덮어쓰기)
# =========================================================
REF_DIR          = os.environ.get("REF_DIR", "NCBI-RefSeq-16s-202505")

# QIIME2 분석 산출 루트
ANALYSIS         = os.environ.get("ANALYSIS", "GutBiomeTech/vsearch_results")

# FASTQ 위치 & 파일 접미사
FASTQ_DIR        = os.environ.get("FASTQ_DIR", "/home/ljb/qiime2_analysis/GutBiomeTech/M17_fastq")
SUFFIX_R1        = os.environ.get("SUFFIX_R1", "_L001_R1_001.fastq.gz")
SUFFIX_R2        = os.environ.get("SUFFIX_R2", "_L001_R2_001.fastq.gz")

# 쓰레드
THREADS_CUTADAPT = os.environ.get("THREADS_CUTADAPT", "24")
THREADS_BLAST    = os.environ.get("THREADS_BLAST", "48")

# MicrobiomeAnalyst 경로/출력
MA_DIR       = os.environ.get("MA_DIR", "GutBiomeTech/MicrobiomeAnalyst")
MA_UPLOAD_DIR = f"{MA_DIR}/upload"  # feature/tax/meta 업로드 파일 위치
MA_OUTDIR    = os.environ.get("MA_OUTDIR", f"{MA_DIR}/results_R_test_MA2")
MA_GROUP     = os.environ.get("MA_GROUP", "Group1")    # 메타데이터의 그룹 열명
RAREFY_Q     = os.environ.get("RAREFY_Q", "0.00")      # 0.00이면 min, 0~1이면 quantile

# LEfSe 전용 출력 디렉토리
LEFSE_OUTDIR = os.environ.get(
    "LEFSE_OUTDIR",
    "/home/ljb/qiime2_analysis/GutBiomeTech/MicrobiomeAnalyst/results_R_test_MA2/04_LEfse",
)

# Picrust2 출력 디렉토리

PICRU    = f"GutBiomeTech/vsearch_results/picrust2"
TGT      = f"{PICRU}/picrust2_data"

# Rscript 바이너리
RS_BIN = os.environ.get("RS_BIN", "/usr/bin/Rscript")

# 내부 파일(센티널 등)
MANIFEST     = f"{ANALYSIS}/sample_manifest.txt"
ANALYZE_DONE = ".done/_3_analyze_sh.done"


REPORT_DIR = "report"

# =========================================================
# 최종 타깃
# =========================================================
rule all:
    input:
        # QIIME2 분석 완료 센티널
        ANALYZE_DONE,
        # MicrobiomeAnalyst 업로드 파일
        f"{MA_UPLOAD_DIR}/feature_table_for_MA.txt",
        f"{MA_UPLOAD_DIR}/metadata_for_MA.txt",
        f"{MA_UPLOAD_DIR}/taxonomy_for_MA.txt",
        # MA 전체 파이프라인 완료
        f"{MA_OUTDIR}/.done_pipeline",
        # pairwise LEfSe 완료
        f"{LEFSE_OUTDIR}/.done_pairwise",
        # 그림 수집 완료
        f"{MA_OUTDIR}/.done_figures",
	# picrust2
        f"{ANALYSIS}/picrust2/.done_picrust2",
        f"{TGT}/KO_unstrat.tsv",
        f"{TGT}/pathway_unstrat.tsv",
        f"{PICRU}/plots_by_class/errorbar_by_class.pdf",
	# Report
        f"{REPORT_DIR}/index.html"

# =========================================================
# 1) RefSeq 16S DB 구축 (이미 있으면 스킵)
# =========================================================
rule build_db:
    output:
        db  = f"{REF_DIR}/ncbi-refseqs-blastdb.qza",
        tax = f"{REF_DIR}/ncbi-refseqs-taxonomy-derep.qza"
    params:
        script  = "rules/2_NCBI_build_db.sh",
        ref_dir = REF_DIR
    shell:
        r"""
        set -euo pipefail
        REF_DIR="{params.ref_dir}" bash "{params.script}"
        """

# =========================================================
# 1.5) FASTQ → sample_manifest.txt 생성
# =========================================================
rule make_manifest:
    input:
        # 순서 보장을 위해 build_db 출력에 의존 (파일 자체를 쓰진 않음)
        db  = rules.build_db.output.db,
        tax = rules.build_db.output.tax
    output:
        MANIFEST
    params:
        script    = "rules/make_manifest.sh",
        analysis  = ANALYSIS,
        fastq_dir = FASTQ_DIR,
        r1        = SUFFIX_R1,
        r2        = SUFFIX_R2
    shell:
        r"""
        set -euo pipefail
        ANALYSIS="{params.analysis}" FASTQ_DIR="{params.fastq_dir}" \
        SUFFIX_R1="{params.r1}" SUFFIX_R2="{params.r2}" \
        MANIFEST="{output}" bash "{params.script}"
        """

# =========================================================
# 2) 샘플 분석 (cutadapt/BLAST 등) — manifest & DB 입력 요구
# =========================================================
rule run_qiime2:
    input:
        manifest = MANIFEST,
        db       = rules.build_db.output.db,
        tax      = rules.build_db.output.tax
    output:
        done_file = ".done/_3_analyze_sh.done",
	ft   = f"{ANALYSIS}/ASV_quantified/feature-table.txt",
        rep_qza   = f"{ANALYSIS}/rep-seqs-dada2-2.qza",
        tab_qza   = f"{ANALYSIS}/table-dada2-2.qza"
    params:
        script        = "rules/3_analyze.sh",
        analysis      = ANALYSIS,
        threads_cut   = THREADS_CUTADAPT,
        threads_blast = THREADS_BLAST,
        ref_dir       = REF_DIR,
        manifest      = MANIFEST
    conda:
        "./envs/qiime2-amplicon-ubuntu-latest-conda.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p .done
        ANALYSIS="{params.analysis}" \
        THREADS_CUTADAPT="{params.threads_cut}" THREADS_BLAST="{params.threads_blast}" \
        REF_DIR="{params.ref_dir}" MANIFEST="{params.manifest}" \
        bash "{params.script}"

        # 산출물 확인
        test -s "{output.ft}"
        test -s "{output.rep_qza}"
        test -s "{output.tab_qza}"

        date > {output.done_file}
        """

# =========================================================
# 3) MicrobiomeAnalyst 업로드용 파일(feature/tax/meta) 생성
# =========================================================
rule make_MicrobiomeAnalystR_input:
    input:
        ANALYZE_DONE,  # 순서 보장용
        f"{ANALYSIS}/ASV_quantified/feature-table.txt",
        f"{ANALYSIS}/sample_manifest.txt",
        script = "rules/make_MicroAnalysis_input.sh"
    output:
        ft   = f"{MA_UPLOAD_DIR}/feature_table_for_MA.txt",
        meta = f"{MA_UPLOAD_DIR}/metadata_for_MA.txt",
        tax  = f"{MA_UPLOAD_DIR}/taxonomy_for_MA.txt"
    params:
        analysis = ANALYSIS,
        ma_dir   = MA_DIR
    shell:
        r"""
        set -euo pipefail
        chmod +x "{input.script}"
        ANALYSIS="{params.analysis}" MA_DIR="{params.ma_dir}" \
        bash "{input.script}"

        # 산출물 존재 확인
        test -s "{output.ft}"
        test -s "{output.meta}"
        test -s "{output.tax}"
        """

# =========================================================
# 4) MicrobiomeAnalystR 전체 파이프라인 (alpha/beta/LEfSe 등)
# Rscript 내에서 library 경로를 지정해야 함
# =========================================================
rule run_MicrobiomeAnalyst_pipeline:
    input:
        ft    = f"{MA_UPLOAD_DIR}/feature_table_for_MA.txt",
        meta  = f"{MA_UPLOAD_DIR}/metadata_for_MA.txt",
        tax   = f"{MA_UPLOAD_DIR}/taxonomy_for_MA.txt",
        script = "scripts/microbiomeanalyst_pipeline.R"
    output:
        touch(f"{MA_OUTDIR}/.done_pipeline")
    params:
        outdir = MA_OUTDIR,
        group  = MA_GROUP,
        rq     = RAREFY_Q,
	pcut   = 0.1,
	count_cut  = 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        "{RS_BIN}" "{input.script}" \
          --feat_fp "{input.ft}" \
          --tax_fp  "{input.tax}" \
          --meta_fp "{input.meta}" \
          --outdir  "{params.outdir}" \
          --group   "{params.group}" \
	  --p_cut   "{params.pcut}" \
          --count_cut   "{params.count_cut}" \
          --rarefy_q {params.rq}
        date > "{output}"
        """

# =========================================================
# 5) pairwise LEfSe (ALL + 그룹쌍) — 전용 폴더로 출력
# Rscript 내에서 library 경로를 지정해야 함

# =========================================================
rule run_MicrobiomeAnalyst_LefSe:
    input:
        ft     = f"{MA_UPLOAD_DIR}/feature_table_for_MA.txt",
        meta   = f"{MA_UPLOAD_DIR}/metadata_for_MA.txt",
        tax    = f"{MA_UPLOAD_DIR}/taxonomy_for_MA.txt",
        script = "scripts/run_ma_lefse.R"
    output:
        touch(f"{LEFSE_OUTDIR}/.done_pairwise")
    params:
        outdir = LEFSE_OUTDIR,
        group  = MA_GROUP,
        rq     = RAREFY_Q,
        pcut   = 0.1,
        count_cut  = 4

    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        "{RS_BIN}" "{input.script}" \
          --feat_fp "$(realpath "{input.ft}")" \
          --tax_fp  "$(realpath "{input.tax}")" \
          --meta_fp "$(realpath "{input.meta}")" \
          --outdir  "$(realpath -m "{params.outdir}")" \
          --group   "{params.group}" \
          --p_cut   "{params.pcut}" \
          --count_cut   "{params.count_cut}" \
          --rarefy_q {params.rq}
        date > "{output}"
        """

# =========================================================
# 6) 그림/결과물 수집 (평탄화)
# =========================================================
rule collect_ma_figures:
    input:
        f"{MA_OUTDIR}/.done_pipeline"
    output:
        touch(f"{MA_OUTDIR}/.done_figures")
    params:
        outdir = MA_OUTDIR,
        dest   = f"{MA_OUTDIR}/99_Figures"
    shell:
        r"""
        set -euo pipefail
        OUT="{params.outdir}"
        DEST="{params.dest}"
        mkdir -p "$DEST"

        # PNG/PDF 수집 (자기 자신 폴더 제외)
        find "$OUT" -type f \( -iname '*.png' -o -iname '*.pdf' \) ! -path "$DEST/*" -print0 \
          | while IFS= read -r -d '' f; do
                rel="${{f#"$OUT/"}}"
                safe="${{rel//\//__}}"
                cp -f "$f" "$DEST/$safe"
            done

        # 인덱스 작성
        printf "collected\toriginal_path\n" > "$DEST/_index.tsv"
        find "$DEST" -type f \( -iname '*.png' -o -iname '*.pdf' \) -printf "%f\t$OUT/99_Figures/%f\n" >> "$DEST/_index.tsv"

        date > "{output}"
        """
# =========================================================
# 7) PICRUSt2: QIIME2 산출물 → PICRUSt2 입력(export) → 실행
#    rep-seqs-dada2-2.qza, table-dada2-2.qza 를 사용
# =========================================================

rule export_repseqs_fna:
    input:
        rep = f"{ANALYSIS}/rep-seqs-dada2-2.qza"
    output:
        fna = f"{ANALYSIS}/ASV_quantified/rep_seqs.fna"
    shell:
        r"""
        set -euo pipefail
        tmpdir="$(mktemp -d)"
        qiime tools export --input-path "{input.rep}" --output-path "$tmpdir"
        mkdir -p "$(dirname "{output.fna}")"
        cp "$tmpdir/dna-sequences.fasta" "{output.fna}"
        """

rule export_table_biom:
    input:
        tab = f"{ANALYSIS}/table-dada2-2.qza"
    output:
        biom = f"{ANALYSIS}/ASV_quantified/table_seqs.biom"
    shell:
        r"""
        set -euo pipefail
        tmpdir="$(mktemp -d)"
        qiime tools export --input-path "{input.tab}" --output-path "$tmpdir"
        mkdir -p "$(dirname "{output.biom}")"
        cp "$tmpdir/feature-table.biom" "{output.biom}"
        """

rule run_picrust2:
    input:
        fna  = f"{ANALYSIS}/ASV_quantified/rep_seqs.fna",
        biom = f"{ANALYSIS}/ASV_quantified/table_seqs.biom"
    output:
        outdir = directory(f"{ANALYSIS}/picrust2"),
        done   = f"{ANALYSIS}/picrust2/.done_picrust2"
    threads: 40
    conda: "envs/picrust2.yml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR="{output.outdir}"
        WORKDIR="$(mktemp -d)"
        trap 'rm -rf "$WORKDIR"' EXIT

        # 임시 작업 디렉토리로 먼저 실행 (기존 OUTDIR가 있어도 영향 없음)
        picrust2_pipeline.py \
          -s "{input.fna}" \
          -i "{input.biom}" \
          -o "$WORKDIR/out" \
          -p {threads}

        # 이전 결과가 있으면 삭제 후 교체
        rm -rf "$OUTDIR"
        mkdir -p "$(dirname "$OUTDIR")"
        mv "$WORKDIR/out" "$OUTDIR"

        date > "{output.done}"
        """

# Bash로 QC + export (unstrat only)
rule qc_and_export_picrust2:
    input:
        flag = f"{ANALYSIS}/picrust2/.done_picrust2"     # PICRUSt2 완료 플래그
    output:
        ko_un = f"{TGT}/KO_unstrat.tsv",
        # EC는 없을 수도 있으니 명시적 산출물에서 제외(원하면 추가)
        pw_un = f"{TGT}/pathway_unstrat.tsv",
        nsti  = f"{TGT}/nsti_summary.txt"
    params:
        OUTDIR = f"{PICRU}",
        TGT    = f"{TGT}"
    shell:
        r"""
        bash scripts/picrust2_qc_and_export.sh "{params.OUTDIR}" "{params.TGT}"
        """


rule plot_kegg_by_class:
    input:
        ko_un = f"{TGT}/KO_unstrat.tsv",
        meta  = f"{MA_DIR}/upload/metadata_for_MA.txt"
    output:
        pdf = f"{PICRU}/plots_by_class/errorbar_by_class.pdf"
    params:
        out_dir = f"{PICRU}/plots_by_class",
        group   = "Group1",
        topn    = 5
    shell:
        r"""
        mkdir -p "{params.out_dir}"
        Rscript scripts/plot_kegg_by_class.R \
          "{input.ko_un}" "{input.meta}" "{params.group}" "{params.out_dir}" "{params.topn}"
        # plot_kegg_by_class.R가 out_dir/errorbar_by_class.pdf를 저장하도록 이미 되어 있어야 함
        test -s "{output.pdf}"
        """

rule build_html_report:
    input:
        ma_done = f"{MA_OUTDIR}/.done_figures",
        lefse_done = f"{LEFSE_OUTDIR}/.done_pairwise",
        pic_done = f"{PICRU}/.done_picrust2",
        # PICRUSt2 플롯 + TSV들 (있으면 모두 의존)
        pic_pdf  = f"{PICRU}/plots_by_class/errorbar_by_class.pdf",
        pic_src  = f"{PICRU}/plots_by_class/errorbar_by_class_source.tsv",
        pic_daa  = f"{PICRU}/plots_by_class/daa_kegg_raw.tsv",
        pic_kw   = f"{PICRU}/plots_by_class/daa_kegg_kruskal.tsv",
        pic_ann  = f"{PICRU}/plots_by_class/daa_kegg_annotated.tsv",
        pic_top  = f"{PICRU}/plots_by_class/errorbar_by_class_top_terms.tsv",
    output:
        html = f"{REPORT_DIR}/index.html"
    params:
        ma_outdir        = MA_OUTDIR,
        ma_fig_dir       = f"{MA_OUTDIR}/99_Figures",
        lefse_csv        = f"{LEFSE_OUTDIR}/lefse_de_output.csv",
        lefse_dir        = LEFSE_OUTDIR,  # overall/pairwise 구조 스캔
        picrust_plot_dir = f"{PICRU}/plots_by_class",
        outdir           = REPORT_DIR,
        py               = SCRIPT_BUILD_REPORT
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"
        python3 "{params.py}" \
          --ma_outdir "{params.ma_outdir}" \
          --ma_fig_dir "{params.ma_fig_dir}" \
          --lefse_csv "{params.lefse_csv}" \
          --lefse_dir "{params.lefse_dir}" \
          --picrust_plot_dir "{params.picrust_plot_dir}" \
          --outdir "{params.outdir}" \
          --mode copy
        test -s "{output.html}"
        """

