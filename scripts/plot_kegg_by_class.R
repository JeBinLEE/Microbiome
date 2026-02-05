#!/usr/bin/env Rscript
.libPaths("/home/ljb/R/x86_64-pc-linux-gnu-library/4.4")

#.libPaths("/home/mysns/miniconda3/envs/microbiome-ma/lib/R/library")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(stringr)
  library(ggplot2); library(ggprism); library(ggpicrust2) ;library(rlang)
  
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_kegg_by_class.R <KO_unstrat.tsv> <metadata.tsv> <group_col> <out_dir> [top_n_per_class=5]")
}
ko_unstrat_tsv <- args[1]
metadata_path  <- args[2]
group_col      <- args[3]
out_dir        <- args[4]

# ko_unstrat_tsv <- "~/microbiome/Microbiome/GutBiomeTech/vsearch_results/picrust2/picrust2_data/KO_unstrat.tsv"
# metadata_path  <- "~/microbiome/Microbiome/GutBiomeTech/MicrobiomeAnalyst/upload/metadata_for_MA.txt"
# group_col      <- "Group1"
# out_dir        <- args[4]

top_n_per_class <- ifelse(length(args) >= 5, as.integer(args[5]), 5)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[picrust] KO_unstrat: ", ko_unstrat_tsv)
message("[picrust] metadata  : ", metadata_path)
message("[picrust] group_col : ", group_col)
message("[picrust] out_dir   : ", out_dir)

## ---------- utils ----------
make_placeholder_pdf <- function(fp, msg = "No pathways to display") {
  # Use base pdf device (available on almost all servers); avoid cairo_pdf dependency.
  grDevices::pdf(fp, width = 7, height = 3)
  graphics::plot.new()
  graphics::text(0.5, 0.5, msg, cex = 1.2)
  grDevices::dev.off()
}

sanitize_ko_unstrat <- function(path_in) {
  # 안전하게 문자로 읽고(주석 '#' 무시), 우리가 타입을 강제 변환
  df <- readr::read_tsv(
    file = path_in,
    comment = "#",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
  
  # 첫 컬럼명 정규화: function / function. -> function.
  if (ncol(df) == 0) stop("Empty KO file: ", path_in)
  if (!names(df)[1] %in% c("function", "function.")) {
    # 드문 케이스: 첫 셀에 "function"이 데이터로 흘러들면 헤더를 재구성
    names(df)[1] <- "function."
  } else {
    names(df)[1] <- "function."
  }
  
  # KO ID 정규화 (대문자 'K'로)
  df[["function."]] <- toupper(df[["function."]])
  
  # KO ID 정합성: K##### 만 허용
  ok_row <- grepl("^K\\d{5}$", df[["function."]])
  dropped <- sum(!ok_row)
  if (dropped > 0) {
    message("[sanitize] Drop ", dropped, " non-KO rows (incl. stray 'function').")
    df <- df[ok_row, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("No valid KO IDs remain after sanitization.")
  
  # 샘플 열을 숫자로 강제 (문자/빈칸/NA → NA_real_ → 0)
  if (ncol(df) > 1) {
    for (j in 2:ncol(df)) {
      v <- df[[j]]
      # 공백/빈문자 처리
      v[trimws(v) %in% c("", "NA", "NaN", "nan", "NAN")] <- NA_character_
      vnum <- suppressWarnings(as.numeric(v))
      vnum[is.na(vnum)] <- 0
      df[[j]] <- vnum
    }
  }
  
  # 임시 파일로 기록 (첫 컬럼이 function. 이어야 함)
  tmp <- tempfile(fileext = ".tsv")
  readr::write_tsv(df, tmp)
  return(tmp)
}

# -------------------------------
# main wrapper: never fail snakemake
# -------------------------------
main <- function() {

## ---------- 1) 입력 로드 ----------
metadata <- read_delim(metadata_path, delim = "\t", show_col_types = FALSE)
stopifnot("Sample" %in% names(metadata))
stopifnot(group_col %in% names(metadata))

# KO 파일 정화 → ggpicrust2로 변환
ko_fixed <- sanitize_ko_unstrat(ko_unstrat_tsv)
kegg_abundance <- ko2kegg_abundance(ko_fixed)

common <- intersect(colnames(kegg_abundance), metadata$Sample)
if (length(common) < 2) stop("Not enough common samples after intersection.")
kegg_abundance <- kegg_abundance[, common, drop = FALSE]
metadata <- metadata |>
  dplyr::filter(Sample %in% common) |>
  dplyr::slice(match(colnames(kegg_abundance), Sample))

## ---------- 2) DAA & annotation ----------
## ---------- 2) DAA & annotation ----------
daa_kegg <- tryCatch({
pathway_daa(
  kegg_abundance,
  metadata,
  group      = group_col,
  daa_method = "ALDEx2"
)
}, error = function(e) {
  msg <- paste0('pathway_daa failed: ', conditionMessage(e))
  message(msg)
  writeLines(msg, file.path(out_dir, 'pathway_daa_ERROR.txt'))
  make_placeholder_pdf(file.path(out_dir, 'errorbar_all.pdf'), msg)
  make_placeholder_pdf(file.path(out_dir, 'errorbar_by_class.pdf'), msg)
  quit(save = 'no', status = 0)
})

daa_kegg_kw <- daa_kegg |>
  dplyr::filter(grepl("Kruskal|KW|kruskal|ALDEx2", .data$method, ignore.case = TRUE))

# If multiple methods remain (e.g. ALDEx2 returns Welch and Wilcoxon), pick the first one
if (nrow(daa_kegg_kw) > 0) {
  methods_found <- unique(daa_kegg_kw$method)
  if (length(methods_found) > 1) {
    chosen <- methods_found[1]
    message("[picrust] Multiple methods found: ", paste(methods_found, collapse=", "), ". Using: ", chosen)
    daa_kegg_kw <- daa_kegg_kw |> dplyr::filter(method == chosen)
  }
}

## (1) DAA 결과가 아예 없으면 → 더미 pdf 만들고 정상 종료
if (is.null(daa_kegg_kw) || nrow(daa_kegg_kw) == 0) {
  message("No KEGG pathways from pathway_daa (ALDEx2/Kruskal). Skipping KEGG class plot.")
  
  # 디버깅용 원본도 남겨두기
  readr::write_tsv(daa_kegg,    file.path(out_dir, "daa_kegg_raw.tsv"))
  readr::write_tsv(daa_kegg_kw, file.path(out_dir, "daa_kegg_kruskal.tsv"))
  
  make_placeholder_pdf(
    file.path(out_dir, "errorbar_all.pdf"),
    "No KEGG pathways from ALDEx2/Kruskal"
  )
  make_placeholder_pdf(
    file.path(out_dir, "errorbar_by_class.pdf"),
    "No KEGG pathways from ALDEx2/Kruskal"
  )
  
  quit(save = "no", status = 0)
}

## (2) KEGG annotation
daa_kegg_annot <- pathway_annotation(
  pathway        = "KO",
  daa_results_df = daa_kegg_kw,
  ko_to_kegg     = TRUE
)

## annotation= 0 → dummi  pdf create
if (is.null(daa_kegg_annot) || nrow(daa_kegg_annot) == 0) {
  message("KEGG annotation returned 0 rows. Skipping KEGG class plot.")
  
  readr::write_tsv(daa_kegg,    file.path(out_dir, "daa_kegg_raw.tsv"))
  readr::write_tsv(daa_kegg_kw, file.path(out_dir, "daa_kegg_kruskal.tsv"))
  
  make_placeholder_pdf(
    file.path(out_dir, "errorbar_all.pdf"),
    "No KEGG annotation rows"
  )
  make_placeholder_pdf(
    file.path(out_dir, "errorbar_by_class.pdf"),
    "No KEGG annotation rows"
  )
  
  quit(save = "no", status = 0)
}

## (3) pathway_class 처리
if (!("pathway_class" %in% names(daa_kegg_annot))) {
  daa_kegg_annot$pathway_class <- NA_character_
}

naa <- is.na(daa_kegg_annot$pathway_class)
if (all(naa)) {
  message("All pathways have NA class. Setting pathway_class = 'Unknown'.")
  daa_kegg_annot$pathway_class <- "Unknown"
} else if (any(naa)) {
  message("[ggpicrust2] drop ", sum(naa), " rows with NA class.")
  daa_kegg_annot <- daa_kegg_annot[!naa, , drop = FALSE]
}

## (4) 여기까지 왔으면 최소 1행 보장 → effect_size 추가해도 안전
eff_candidates <- c(
  "effect", "diff.btw", "diff.between",
  "logFC", "cohen_d", "AUC", "stat", "W"
)
hit <- eff_candidates[eff_candidates %in% names(daa_kegg_annot)]
daa_kegg_annot$effect_size <- NA_real_
if (length(hit)) {
  daa_kegg_annot$effect_size <- suppressWarnings(
    as.numeric(daa_kegg_annot[[hit[1]]])
  )
}

## ---------- 2.1) save source data----------
readr::write_tsv(
  tibble::as_tibble(
    tibble::rownames_to_column(as.data.frame(kegg_abundance), "feature")
  ),
  file.path(out_dir, "errorbar_by_class_source.tsv")
)
readr::write_tsv(daa_kegg,       file.path(out_dir, "daa_kegg_raw.tsv"))
readr::write_tsv(daa_kegg_kw,    file.path(out_dir, "daa_kegg_kruskal.tsv"))
readr::write_tsv(daa_kegg_annot, file.path(out_dir, "daa_kegg_annotated.tsv"))

top_terms <- daa_kegg_annot |>
  dplyr::group_by(pathway_class) |>
  dplyr::arrange(p_adjust, .by_group = TRUE) |>
  dplyr::slice_head(n = top_n_per_class) |>
  dplyr::ungroup() |>
  dplyr::select(pathway_class, feature, pathway_name, p_adjust, effect_size)

readr::write_tsv(
  top_terms,
  file.path(out_dir, "errorbar_by_class_top_terms.tsv")
)

## ---------- 3) 플로팅 ----------
has_rows <- nrow(daa_kegg_annot) > 0

# 전체
if (has_rows) {
  p_all <- pathway_errorbar(
    abundance          = kegg_abundance,
    daa_results_df     = daa_kegg_annot,
    Group              = metadata[[group_col]],
    p_values_threshold = 0.05,
    order              = "pathway_class",
    select             = NULL,
    ko_to_kegg         = TRUE,
    p_value_bar        = FALSE,
    colors             = NULL,
    x_lab              = "pathway_name"
  )
  ggsave(file.path(out_dir, "errorbar_all.pdf"), p_all,
         width = 12, height = 15, units = "in", device = "pdf")
} else {
  make_placeholder_pdf(file.path(out_dir, "errorbar_all.pdf"),
                       "No significant pathways (global)")
}

# class별
plot_list <- list()
if (has_rows) {
  classes <- daa_kegg_annot |>
    dplyr::filter(!is.na(pathway_class)) |>
    dplyr::pull(pathway_class) |>
    unique()
  
  for (cls in classes) {
    df_cls <- dplyr::filter(daa_kegg_annot, pathway_class == cls)
    sel_features <- if (is.finite(top_n_per_class)) {
      df_cls |> dplyr::arrange(p_adjust) |> dplyr::slice_head(n = top_n_per_class) |> dplyr::pull(feature)
    } else df_cls$feature
    if (!length(sel_features)) next
    df_cls <- dplyr::filter(df_cls, feature %in% sel_features)
    
    p <- pathway_errorbar(
      abundance          = kegg_abundance,
      daa_results_df     = df_cls,
      Group              = metadata[[group_col]],
      p_values_threshold = 0.05,
      order              = "group",
      select             = sel_features,
      ko_to_kegg         = TRUE,
      p_value_bar        = FALSE,
      colors             = NULL,
      x_lab              = "pathway_name"
    )
    
    plot_list[[cls]] <- p
    ggsave(file.path(out_dir, paste0("errorbar_", gsub("[^A-Za-z0-9]+","_", cls), ".pdf")),
           p, width = 12, height = max(6, 0.42 * length(sel_features)),
           units = "in", device = "pdf")
  }
}

if (length(plot_list)) {
  pdf(file.path(out_dir, "errorbar_by_class.pdf"), width = 12, height = 7)
  invisible(lapply(plot_list, print))
  dev.off()
} else {
  make_placeholder_pdf(file.path(out_dir, "errorbar_by_class.pdf"),
                       "No class-level plots")
}

message("[OK] PICRUSt2 KEGG plots exported to: ", out_dir)
}

tryCatch(
  main(),
  error = function(e) {
    msg <- paste0("[ERROR] plot_kegg_by_class.R failed: ", conditionMessage(e))
    message(msg)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(msg, file.path(out_dir, "PLOT_KEGG_ERROR.txt"))
    make_placeholder_pdf(file.path(out_dir, "errorbar_all.pdf"), msg)
    make_placeholder_pdf(file.path(out_dir, "errorbar_by_class.pdf"), msg)
    quit(save = "no", status = 0)
  }
)

