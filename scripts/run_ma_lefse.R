#!/usr/bin/env Rscript
.libPaths("/home/ljb/R/x86_64-pc-linux-gnu-library/4.4")

suppressPackageStartupMessages({
  library(optparse)
  library(MicrobiomeAnalystR)
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(rlang)
})

## -------------------------------
## CLI 옵션
## -------------------------------
option_list <- list(
  make_option("--feat_fp",  type="character", help="feature_table_for_MA.txt"),
  make_option("--tax_fp",   type="character", help="taxonomy_for_MA.txt"),
  make_option("--meta_fp",  type="character", help="metadata_for_MA.txt"),
  make_option("--outdir",   type="character", help="Output directory"),
  make_option("--group",    type="character", default="Group1"),
  make_option("--p_cut",    type="double",    default=0.10),
  make_option("--count_cut",type="double",    default=4),
  make_option("--rarefy_q", type="double",    default=0.00,
              help="0.00=min depth; (0,1]=quantile for rarefaction depth (default: %default)")
)
opt <- parse_args(OptionParser(option_list = option_list))
req <- c("feat_fp","tax_fp","meta_fp","outdir")
if (any(vapply(opt[req], is.null, TRUE))) {
  stop("Missing required args. Use: --feat_fp --tax_fp --meta_fp --outdir [--group] [--rarefy_q]")
}

## 절대경로
feat_fp   <- normalizePath(opt$feat_fp, mustWork=TRUE)
tax_fp    <- normalizePath(opt$tax_fp,  mustWork=TRUE)
meta_fp   <- normalizePath(opt$meta_fp, mustWork=TRUE)
outdir    <- normalizePath(opt$outdir,  mustWork=FALSE)
group_col <- opt$group
rarefy_q  <- opt$rarefy_q
p_cut     <- opt$p_cut
count_cut <- opt$count_cut

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
setwd(outdir)

message("== MicrobiomeAnalystR end-to-end + pairwise LEfSe (Genus) ==")
message("Output dir: ", outdir)
message("Group column: ", group_col, " | rarefy_q: ", rarefy_q)


# -------------------------------
# utils: never fail the pipeline
# -------------------------------
make_placeholder <- function(outdir, msg) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  # human-readable reason
  writeLines(as.character(msg), file.path(outdir, "LEFSE_SKIPPED.txt"))
  # also create an empty CSV so downstream tooling won't crash
  csv_fp <- file.path(outdir, "lefse_de_output.csv")
  if (!file.exists(csv_fp)) {
    writeLines("feature,lda,p_value,p_adjust,group", csv_fp)
  }
}

safe_quit_ok <- function(outdir, msg) {
  message(msg)
  try(make_placeholder(outdir, msg), silent = TRUE)
  quit(save = "no", status = 0)
}



# -------------------------------
# 0) Quick metadata sanity check (LEfSe needs >=2 groups and >=3 replicates per group)
# -------------------------------
meta_chk <- tryCatch(readr::read_tsv(meta_fp, show_col_types = FALSE), error=function(e) NULL)
if (is.null(meta_chk)) {
  msg <- paste0('Skipping LEfSe: failed to read metadata: ', meta_fp)
  message(msg)
  make_placeholder(outdir, msg)
  quit(save='no', status=0)
}
if (!(group_col %in% names(meta_chk))) {
  msg <- paste0('Skipping LEfSe: group column not found in metadata: ', group_col)
  message(msg)
  make_placeholder(outdir, msg)
  quit(save='no', status=0)
}
# Some MA metadata use first column as Sample ID
if (!("Sample" %in% names(meta_chk))) {
  names(meta_chk)[1] <- "Sample"
}
meta_chk[[group_col]] <- as.character(meta_chk[[group_col]])
meta_chk <- meta_chk[!is.na(meta_chk[[group_col]]) & meta_chk[[group_col]] != "", , drop=FALSE]
ct <- table(meta_chk[[group_col]])
if (length(ct) < 2) {
  msg <- paste0('Skipping LEfSe: need >=2 groups. Found: ', paste(names(ct), collapse=','))
  message(msg)
  make_placeholder(outdir, msg)
  quit(save='no', status=0)
}
if (any(ct < 3)) {
  msg <- paste0('Skipping LEfSe: each group needs >=3 replicates. Counts: ', paste(names(ct), ct, sep='=', collapse=', '))
  message(msg)
  make_placeholder(outdir, msg)
  quit(save='no', status=0)
}
## 유틸: 파일명 안전 토큰
safe_token <- function(x) gsub("[^A-Za-z0-9]+","_", x)

# -------------------------------
# condenseOTUs 자동 탐색/로드 (없으면 shim)
# -------------------------------
find_and_load_condense <- function(
    roots = c(
      getwd(),
      dirname(getwd()),
      "/home/ljb/qiime2_analysis/GutBiomeTech/MicrobiomeAnalyst",
      "/home/ljb/qiime2_analysis/GutBiomeTech/MicrobiomeAnalyst/MicrobiomeAnalystR",
      system.file(package = "MicrobiomeAnalystR")
    )
){
  roots <- unique(normalizePath(roots[file.exists(roots)], mustWork = FALSE))
  hits  <- character(0)
  for (rt in roots) {
    cand_dirs <- c(file.path(rt, "rscripts"), rt)
    cand_dirs <- cand_dirs[file.exists(cand_dirs)]
    for (dd in cand_dirs) {
      rfiles <- list.files(dd, pattern = "[.]R$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      if (!length(rfiles)) next
      hit <- rfiles[ vapply(rfiles, function(f){
        any(grepl("\\bcondenseOTUs\\s*<-\\s*function", readLines(f, warn = FALSE), perl = TRUE))
      }, logical(1)) ]
      if (length(hit)) {
        src_dir <- dirname(hit[1])
        message("Loading rscripts from: ", src_dir)
        rs <- list.files(src_dir, pattern = "[.]R$", full.names = TRUE, recursive = FALSE)
        invisible(lapply(rs, function(f){ message("  source: ", basename(f)); source(f, chdir = TRUE) }))
        hits <- c(hits, hit[1]); break
      }
    }
    if (length(hits)) break
  }
  exists("condenseOTUs", mode = "function")
}
ok <- find_and_load_condense()
message("condenseOTUs found? ", ok)
if (!ok) {
  message("Injecting a minimal shim for condenseOTUs() to unblock normalization…")
  condenseOTUs <- function(otab2, mode = "tax") {
    if (is.data.frame(otab2)) otab2 <- as.matrix(otab2)
    mode(otab2) <- "numeric"
    return(otab2)
  }
}
stopifnot(exists("condenseOTUs", mode = "function"))

# -------------------------------
# 1) 베이스 파이프라인 (전체)
# -------------------------------
mbSet <- Init.mbSetObj()
mbSet <- SetModuleType(mbSet, "mdp")

mbSet <- tryCatch(ReadSampleTable(mbSet, meta_fp), error=function(e){ msg<-paste0('LEfSe: ReadSampleTable failed: ', conditionMessage(e)); message(msg); make_placeholder(outdir,msg); quit(save='no',status=0) })
mbSet <- tryCatch(Read16STaxaTable(mbSet, tax_fp), error=function(e){ msg<-paste0('LEfSe: Read16STaxaTable failed: ', conditionMessage(e)); message(msg); make_placeholder(outdir,msg); quit(save='no',status=0) })
mbSet <- tryCatch(Read16SAbundData(mbSet, feat_fp, "text", "QIIME", "T", "false"), error=function(e){ msg<-paste0('LEfSe: Read16SAbundData failed: ', conditionMessage(e)); message(msg); make_placeholder(outdir,msg); quit(save='no',status=0) })

# Some MA builds return a broken object (atomic vector) without throwing.
if (!is.list(mbSet) || is.null(mbSet$dataSet)) {
  safe_quit_ok(outdir, "Skipping LEfSe: MicrobiomeAnalystR returned an invalid mbSet object after reading inputs")
}

mbSet <- tryCatch(SanityCheckData(mbSet, "text", "sample", "true"),
                  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: SanityCheckData failed: ', conditionMessage(e))) })
mbSet <- tryCatch(SanityCheckSampleData(mbSet),
                  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: SanityCheckSampleData failed: ', conditionMessage(e))) })

# Set group column as the primary metadata.
mbSet <- tryCatch(SetMetaAttributes(mbSet, "1"),
                  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: SetMetaAttributes failed: ', conditionMessage(e))) })

mbSet <- tryCatch(PlotLibSizeView(mbSet, "norm_libsizes_0", "png"),
                  error=function(e) { message('[warn] PlotLibSizeView failed: ', conditionMessage(e)); mbSet })
mbSet <- tryCatch(CreatePhyloseqObj(mbSet, "text", "QIIME", "F", "false"),
                  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: CreatePhyloseqObj failed: ', conditionMessage(e))) })

# -------------------------------
# 2) 전체 필터 → 정규화
# -------------------------------
mbSet <- tryCatch(ApplyAbundanceFilter(mbSet, "prevalence", count_cut, 0.10),
                  error=function(e) { message('[warn] ApplyAbundanceFilter failed: ', conditionMessage(e)); mbSet })
mbSet <- tryCatch(ApplyVarianceFilter(mbSet, "iqr", 0.10),
                  error=function(e) { message('[warn] ApplyVarianceFilter failed: ', conditionMessage(e)); mbSet })

libsizes <- tryCatch(phyloseq::sample_sums(mbSet$dataSet$proc.phyobj), error=function(e) numeric(0))
if (!length(libsizes)) safe_quit_ok(outdir, "Skipping LEfSe: cannot compute library sizes")

depth <- if (rarefy_q > 0 && rarefy_q <= 1) {
  as.integer(stats::quantile(libsizes, probs = rarefy_q, na.rm = TRUE))
} else {
  min(libsizes)
}
depth <- max(1L, min(depth, min(libsizes)))

mbSet <- tryCatch(
  PerformNormalization(
    mbSet,
    "none",     # norm.method
    "colsum",   # scale.method
    "none",     # trans.method
    "true",     # rarefy
    depth
  ),
  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: PerformNormalization failed: ', conditionMessage(e))) }
)

# -------------------------------
# 3) LEfSe (전체 그룹) - Genus
# -------------------------------
mbSet <- tryCatch(
  PerformLefseAnal(mbSet, p_cut, "fdr", 2.0, group_col, "F", "NA", "Genus"),
  error=function(e) { safe_quit_ok(outdir, paste0('Skipping LEfSe: PerformLefseAnal failed: ', conditionMessage(e))) }
)

## ✅ 전체 결과 CSV 저장
rt_all <- mbSet$analSet$lefse$resTable
if (!is.null(rt_all) && NROW(rt_all) > 0) {
  # LDA 숫자 보정
  if ("lda" %in% names(rt_all) && !is.numeric(rt_all$lda)) {
    lda_num <- suppressWarnings(as.numeric(as.character(rt_all$lda)))
    if (any(is.finite(lda_num))) rt_all$lda <- lda_num
  }
  readr::write_csv(rt_all, "lefse_de_output.csv")
} else {
  # 결과가 없으면 빈 CSV라도 생성(파이프라인 호환성용)
  writeLines("feature,lda,p_value,p_adjust,group", "lefse_de_output.csv")
}

## 그림 저장
mbSet <- tryCatch(PlotLEfSeSummary(mbSet, 30, "dot", "lefse_dot_genus", "png"),
                  error=function(e) { message('[warn] PlotLEfSeSummary(dot) failed: ', conditionMessage(e)); mbSet })
mbSet <- tryCatch(PlotLEfSeSummary(mbSet, 30, "bar", "lefse_bar_genus", "png"),
                  error=function(e) { message('[warn] PlotLEfSeSummary(bar) failed: ', conditionMessage(e)); mbSet })

# -------------------------------
# 4) 그룹 페어별 — 없으면 pass
# -------------------------------
meta <- tryCatch(readr::read_table(meta_fp, show_col_types = FALSE), error=function(e) NULL)
if (is.null(meta)) {
  safe_quit_ok(outdir, paste0('Skipping pairwise LEfSe: failed to read metadata: ', meta_fp))
}
# In MA metadata, the first column is usually the sample id (#NAME).
if (!("Sample" %in% names(meta))) names(meta)[1] <- "Sample"
if (!(group_col %in% names(meta))) {
  safe_quit_ok(outdir, paste0('Skipping pairwise LEfSe: group column not found: ', group_col))
}

pairs <- combn(unique(as.character(meta[[group_col]])), 2, simplify = FALSE)

for (pp in pairs) {
  mbSet_filter <- mbSet
  g1 <- pp[1]
  g2 <- pp[2]
  comp_id <- paste0(safe_token(g1), "vs", safe_token(g2))
  message("== Pairwise: ", g1, " vs ", g2, " ==")
  
  # smpl.nm.vec 전역 주입 (UpdateSampleItems가 전역을 읽음)
  smpl.nm.vec <- meta %>% dplyr::filter(.data[[group_col]] %in% c(g1, g2)) %>% dplyr::pull(Sample)
  assign("smpl.nm.vec", smpl.nm.vec, envir = .GlobalEnv)
  
  # 샘플 반영; 실패하면 패스
  mbSet_filter <- tryCatch(UpdateSampleItems(mbSet_filter), error = function(e) { message("[skip] UpdateSampleItems: ", e$message); return(NULL) })
  if (is.null(mbSet_filter)) next
  
  # 샘플이 없으면 패스
  if (is.null(mbSet_filter$dataSet$proc.phyobj) ||
      phyloseq::nsamples(mbSet_filter$dataSet$proc.phyobj) == 0) {
    message("[skip] no samples after filtering: ", g1, "vs", g2); next
  }
  
  # 필터; 실패해도 계속 진행
  mbSet_filter <- tryCatch(ApplyAbundanceFilter(mbSet_filter, "prevalence", count_cut, 0.10),
                           error = function(e) { message("[warn] ApplyAbundanceFilter: ", e$message); return(mbSet_filter) })
  mbSet_filter <- tryCatch(ApplyVarianceFilter(mbSet_filter, "iqr", 0.10),
                           error = function(e) { message("[warn] ApplyVarianceFilter: ", e$message); return(mbSet_filter) })
  
  ## --- Rarefaction 깊이 선택 & 정규화 ---
  libsizes <- tryCatch(phyloseq::sample_sums(mbSet_filter$dataSet$proc.phyobj),
                       error = function(e) numeric(0))
  if (!length(libsizes)) { message("[skip] no libsizes: ", g1, "vs", g2); next }
  min_depth <- min(libsizes)
  
  depth2 <- if (rarefy_q > 0 && rarefy_q <= 1) {
    as.integer(stats::quantile(libsizes, probs = rarefy_q, na.rm = TRUE))
  } else {
    min(libsizes)
  }
  depth2 <- max(1L, min(depth2, min(libsizes)))

  mbSet_filter <- tryCatch(
    PerformNormalization(
      mbSet_filter,
      "none", "colsum", "none",
      "true", depth2
    ),
    error = function(e) { message("[skip] PerformNormalization: ", e$message); return(NULL) }
  )
  if (is.null(mbSet_filter)) next
  
  # LEfSe; 실패하면 패스
  mbSet_filter <- tryCatch(
    PerformLefseAnal(mbSet_filter, p_cut, "fdr", 2.0, group_col, "F", "NA", "Genus"),
    error = function(e) { message("[skip] PerformLefseAnal: ", e$message); return(NULL) }
  )
  if (is.null(mbSet_filter)) next
  
  # 결과 없으면 패스
  rt <- mbSet_filter$analSet$lefse$resTable
  if (is.null(rt) || NROW(rt) == 0) { message("[skip] empty LEfSe result: ", g1, "vs", g2); next }
  
  # LDA 숫자 아닌 경우 유효값만 반영
  if ("lda" %in% names(rt) && !is.numeric(rt$lda)) {
    lda_num <- suppressWarnings(as.numeric(as.character(rt$lda)))
    keep <- is.finite(lda_num) & !is.na(lda_num)
    if (!any(keep)) { message("[skip] non-finite LDA: ", g1, "vs", g2); next }
    mbSet_filter$analSet$lefse$resTable$lda <- lda_num
    rt <- mbSet_filter$analSet$lefse$resTable
  }
  
  ## ✅ 페어별 CSV 저장
  csv_name <- paste0("lefse_", comp_id, ".csv")
  readr::write_csv(rt, csv_name)
  
  # 플롯; 실패해도 건너뜀
  mbSet_filter <- tryCatch(
    PlotLEfSeSummary(mbSet_filter, 30, "dot", paste0("lefse_dot_genus_",g1,"vs",g2), "png"),
    error = function(e) { message("[skip] Plot dot: ", e$message); return(mbSet_filter) }
  )
  mbSet_filter <- tryCatch(
    PlotLEfSeSummary(mbSet_filter, 30, "bar", paste0("lefse_bar_genus_",g1,"vs",g2), "png"),
    error = function(e) { message("[skip] Plot bar: ", e$message); return(mbSet_filter) }
  )
}

message("\n== Done ==\nOutput dir: ", outdir)
