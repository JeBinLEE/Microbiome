#!/usr/bin/env Rscript

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
top_n_per_class <- ifelse(length(args) >= 5, as.integer(args[5]), 5)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[picrust] KO_unstrat: ", ko_unstrat_tsv)
message("[picrust] metadata  : ", metadata_path)
message("[picrust] group_col : ", group_col)
message("[picrust] out_dir   : ", out_dir)

## ---------- utils ----------
make_placeholder_pdf <- function(fp, msg = "No pathways to display") {
  g <- ggplot() + annotate("text", x=0, y=0, label=msg, size=6) +
    xlim(-1,1) + ylim(-1,1) + theme_void()
  ggsave(fp, g, width=6, height=3, device=cairo_pdf)
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
daa_kegg <- pathway_daa(kegg_abundance, metadata, group = group_col, daa_method = "ALDEx2")
daa_kegg_kw <- daa_kegg |> dplyr::filter(grepl("Kruskal|KW|kruskal", .data$method, ignore.case = TRUE))
daa_kegg_annot <- pathway_annotation(pathway = "KO", daa_results_df = daa_kegg_kw, ko_to_kegg = TRUE)

# 유효 class만
naa <- is.na(daa_kegg_annot$pathway_class)
if (any(naa)) {
  message("[ggpicrust2] drop ", sum(naa), " rows with NA class.")
  daa_kegg_annot <- daa_kegg_annot[!naa, , drop = FALSE]
}

# 효과크기 통일: effect_size 생성
eff_candidates <- c("effect","diff.btw","diff.between","logFC","cohen_d","AUC","stat","W")
hit <- eff_candidates[eff_candidates %in% names(daa_kegg_annot)]
daa_kegg_annot$effect_size <- NA_real_
if (length(hit)) {
  daa_kegg_annot$effect_size <- suppressWarnings(as.numeric(daa_kegg_annot[[hit[1]]]))
}

## ---------- 2.1) 소스 데이터 저장 ----------
readr::write_tsv(
  tibble::as_tibble(tibble::rownames_to_column(as.data.frame(kegg_abundance), "feature")),
  file.path(out_dir, "errorbar_by_class_source.tsv")
)
readr::write_tsv(daa_kegg,        file.path(out_dir, "daa_kegg_raw.tsv"))
readr::write_tsv(daa_kegg_kw,     file.path(out_dir, "daa_kegg_kruskal.tsv"))
readr::write_tsv(daa_kegg_annot,  file.path(out_dir, "daa_kegg_annotated.tsv"))

top_terms <- daa_kegg_annot |>
  dplyr::filter(!is.na(pathway_class)) |>
  dplyr::group_by(pathway_class) |>
  dplyr::arrange(p_adjust, .by_group = TRUE) |>
  dplyr::slice_head(n = top_n_per_class) |>
  dplyr::ungroup() |>
  dplyr::select(pathway_class, feature, pathway_name, p_adjust, effect_size)
readr::write_tsv(top_terms, file.path(out_dir, "errorbar_by_class_top_terms.tsv"))

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
         width = 12, height = 15, units = "in", device = cairo_pdf)
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
           units = "in", device = cairo_pdf)
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