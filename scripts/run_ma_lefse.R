#!/usr/bin/env Rscript
.libPaths("/home/ljb/R/x86_64-pc-linux-gnu-library/4.4")
suppressPackageStartupMessages({
  library(optparse)
  library(MicrobiomeAnalystR)
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# -------------------------------
# CLI 옵션 (필수 4개만)
# -------------------------------
option_list <- list(
  make_option("--feat_fp", type="character", help="feature_table_for_MA.txt"),
  make_option("--tax_fp",  type="character", help="taxonomy_for_MA.txt"),
  make_option("--meta_fp", type="character", help="metadata_for_MA.txt"),
  make_option("--outdir",  type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (any(vapply(list(opt$feat_fp, opt$tax_fp, opt$meta_fp, opt$outdir), is.null, TRUE))) {
  stop("Missing required args. Usage: run_ma_lefse.R --feat_fp ... --tax_fp ... --meta_fp ... --outdir ...")
}

# 절대경로 고정 (setwd 전)
feat_fp <- normalizePath(opt$feat_fp, mustWork = TRUE)
tax_fp  <- normalizePath(opt$tax_fp,  mustWork = TRUE)
meta_fp <- normalizePath(opt$meta_fp, mustWork = TRUE)
outdir  <- normalizePath(opt$outdir,  mustWork = FALSE)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

message("== MicrobiomeAnalystR end-to-end + pairwise LEfSe (Genus) ==")
message("Output dir: ", outdir)

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

mbSet <- ReadSampleTable(mbSet, meta_fp)         # Group1 필수
mbSet <- Read16STaxaTable(mbSet, tax_fp)         # Kingdom~Species
mbSet <- Read16SAbundData(mbSet, feat_fp, "text", "QIIME", "T", "false")

mbSet <- SanityCheckData(mbSet, "text", "sample", "true")
mbSet <- SanityCheckSampleData(mbSet)

mbSet <- SetMetaAttributes(mbSet, "1")           # 첫 factor = Group1

mbSet <- PlotLibSizeView(mbSet, "norm_libsizes_0", "png")
mbSet <- CreatePhyloseqObj(mbSet, "text", "QIIME", "F", "false")

# -------------------------------
# 2) 전체 필터 → 정규화
# -------------------------------
mbSet <- ApplyAbundanceFilter(mbSet, "prevalence", 4, 0.20)
mbSet <- ApplyVarianceFilter(mbSet, "iqr", 0.10)

libsizes  <- phyloseq::sample_sums(mbSet$dataSet$proc.phyobj)
min_depth <- min(libsizes)
mbSet <- PerformNormalization(
  mbSet,
  "none",                  # norm.method
  "colsum",                # scale.method
  "none",                  # trans.method
  "true",                  # rarefy
  as.integer(min_depth)    # rarefy.depth
)

# -------------------------------
# 3) LEfSe (전체 그룹) - Genus
# -------------------------------
mbSet <- PerformLefseAnal(mbSet, 0.10, "fdr", 2.0, "Group1", "F", "NA", "Genus")
mbSet <- PlotLEfSeSummary(mbSet, 30, "dot", "lefse_dot_genus", "png")
mbSet <- PlotLEfSeSummary(mbSet, 30, "bar", "lefse_bar_genus", "png")

# -------------------------------
# 4) 그룹 페어별 — 없으면 pass
# -------------------------------
meta <- readr::read_table(meta_fp, show_col_types = FALSE)
if (!("Group1" %in% names(meta))) stop("metadata에 Group1 열이 없습니다.")
if (!("Sample" %in% names(meta))) stop("metadata에 Sample 열이 없습니다.")

pairs <- combn(unique(meta$Group1), 2, simplify = FALSE)

for (pp in pairs) {
  mbSet_filter <- mbSet
  g1 <- pp[1]
  g2 <- pp[2]
  print(paste0(g1,"vs", g2))
  
  # smpl.nm.vec 전역 주입 (UpdateSampleItems가 전역을 읽음)
  smpl.nm.vec <- meta %>% dplyr::filter(Group1 %in% c(g1,g2)) %>% dplyr::pull(Sample)
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
  mbSet_filter <- tryCatch(ApplyAbundanceFilter(mbSet_filter, "prevalence", 4, 0.20),
                           error = function(e) { message("[warn] ApplyAbundanceFilter: ", e$message); return(mbSet_filter) })
  mbSet_filter <- tryCatch(ApplyVarianceFilter(mbSet_filter, "iqr", 0.10),
                           error = function(e) { message("[warn] ApplyVarianceFilter: ", e$message); return(mbSet_filter) })
  
  ## --- Rarefaction 깊이 선택 & 정규화 ---
  libsizes <- tryCatch(phyloseq::sample_sums(mbSet_filter$dataSet$proc.phyobj),
                       error = function(e) numeric(0))
  if (!length(libsizes)) { message("[skip] no libsizes: ", g1, "vs", g2); next }
  min_depth <- min(libsizes)
  
  mbSet_filter <- tryCatch(
    PerformNormalization(
      mbSet_filter,
      "none", "colsum", "none",
      "true", as.integer(min_depth)
    ),
    error = function(e) { message("[skip] PerformNormalization: ", e$message); return(NULL) }
  )
  if (is.null(mbSet_filter)) next
  
  # LEfSe; 실패하면 패스
  mbSet_filter <- tryCatch(
    PerformLefseAnal(mbSet_filter, 0.10, "fdr", 2.0, "Group1", "F", "NA", "Genus"),
    error = function(e) { message("[skip] PerformLefseAnal: ", e$message); return(NULL) }
  )
  if (is.null(mbSet_filter)) next
  
  # 결과 없으면 패스
  rt <- mbSet_filter$analSet$lefse$resTable
  if (is.null(rt) || NROW(rt) == 0) { message("[skip] empty LEfSe result: ", g1, "vs", g2); next }
  
  # LDA 숫자 아닌 경우 유효값 없으면 패스
  if ("lda" %in% names(rt) && !is.numeric(rt$lda)) {
    lda_num <- suppressWarnings(as.numeric(as.character(rt$lda)))
    keep <- is.finite(lda_num) & !is.na(lda_num)
    if (!any(keep)) { message("[skip] non-finite LDA: ", g1, "vs", g2); next }
    mbSet_filter$analSet$lefse$resTable$lda <- lda_num
  }
  
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
