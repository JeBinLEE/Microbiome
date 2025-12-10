#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(withr)
  library(MicrobiomeAnalystR)
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

## ===== CLI =====
opt_list <- list(
  make_option("--feat_fp",  type="character"),
  make_option("--tax_fp",   type="character"),
  make_option("--meta_fp",  type="character"),
  make_option("--outdir",   type="character"),
  make_option("--group",    type="character", default="Group1"),
  make_option("--rank",     type="character", default="Genus"),
  make_option("--lda_cut",  type="double",    default=2.0),
  make_option("--p_cut",    type="double",    default=0.10),
  make_option("--count_cut",    type="double",    default=4),
  make_option("--rarefy_q", type="double",    default=0.00)
)
opt <- parse_args(OptionParser(option_list=opt_list))
stopifnot(!is.null(opt$feat_fp), !is.null(opt$tax_fp), !is.null(opt$meta_fp), !is.null(opt$outdir))

feat_fp <- opt$feat_fp
tax_fp  <- opt$tax_fp
meta_fp <- opt$meta_fp
outdir  <- opt$outdir
grp     <- opt$group
rank_use<- opt$rank
lda_cut <- opt$lda_cut
p_cut   <- opt$p_cut
count_cut   <- opt$count_cut
rq      <- opt$rarefy_q

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
mkdir <- function(name){ d <- file.path(outdir, name); dir.create(d, FALSE, TRUE); d }

## ===== MicrobiomeAnalystR: Load & QC =====
mbSet <- Init.mbSetObj()
mbSet <- SetModuleType(mbSet, "mdp")
mbSet <- ReadSampleTable(mbSet, meta_fp)
mbSet <- Read16STaxaTable(mbSet, tax_fp)
mbSet <- Read16SAbundData(mbSet, feat_fp, "text", "QIIME", "T", "false")
mbSet <- SanityCheckData(mbSet, "text", "sample", "true")
mbSet <- SanityCheckSampleData(mbSet)
mbSet <- SetMetaAttributes(mbSet, "1")

dir_qc <- mkdir("00_QC_Normalization")
with_dir(dir_qc, { try(PlotLibSizeView(mbSet, "norm_libsizes_before", "png"), silent=TRUE) })

mbSet <- CreatePhyloseqObj(mbSet, "text", "QIIME", "F", "false")
mbSet <- ApplyAbundanceFilter(mbSet, "prevalence", count_cut, p_cut) # change to input parameter
mbSet <- ApplyVarianceFilter(mbSet, "iqr", 0.10)

libsizes <- phyloseq::sample_sums(mbSet$dataSet$proc.phyobj)
if (length(libsizes) == 0) stop("No samples detected after filtering.")
min_depth <- if (rq > 0 && rq < 1) floor(quantile(libsizes, rq)) else min(libsizes)

mbSet <- PerformNormalization(mbSet, "none", "colsum", "none", "true", 30000) # change to input parameter
with_dir(dir_qc, { try(PlotLibSizeView(mbSet, "norm_libsizes_after", "png"), silent=TRUE) })

## ===== common ps/meta =====
ps <- mbSet$dataSet$proc.phyobj
if (is.null(ps) || !inherits(ps, "phyloseq")) stop("proc.phyobj is not a phyloseq object.")
meta_df <- as(sample_data(ps), "data.frame")
stopifnot(grp %in% colnames(meta_df))
meta_df[[grp]] <- factor(meta_df[[grp]])

## 작은 헬퍼: 문자열로 formula 생성
fml <- function(response, group_col) reformulate(group_col, response=response)

## ===== 01 Alpha =====
dir_alpha <- mkdir("01_Alpha")
with_dir(dir_alpha, {
  alpha_df <- estimate_richness(ps, measures=c("Chao1","Shannon")) %>%
    tibble::rownames_to_column("Sample") %>%
    left_join(meta_df %>% tibble::rownames_to_column("Sample") %>% select(Sample, all_of(grp)),
              by="Sample")
  write_tsv(alpha_df, "alpha_values.tsv")
  
  ## Kruskal
  kw_chao <- try(kruskal.test(formula = fml("Chao1", grp),   data=alpha_df), silent=TRUE)
  kw_shan <- try(kruskal.test(formula = fml("Shannon", grp), data=alpha_df), silent=TRUE)
  if (!inherits(kw_chao,"try-error")) capture.output(kw_chao, file="kruskal_chao1.txt")
  if (!inherits(kw_shan,"try-error")) capture.output(kw_shan, file="kruskal_shannon.txt")
  
  ## ANOVA (aov)
  aov_chao <- aov(formula = fml("Chao1", grp),   data=alpha_df)
  aov_shan <- aov(formula = fml("Shannon", grp), data=alpha_df)
  sm_chao <- summary(aov_chao)[[1]]; sm_shan <- summary(aov_shan)[[1]]
  lab_chao <- sprintf("P : %.4g; [ANOVA]F-value : %.3f", sm_chao$`Pr(>F)`[1], sm_chao$`F value`[1])
  lab_shan <- sprintf("P : %.4g; [ANOVA]F-value : %.3f", sm_shan$`Pr(>F)`[1], sm_shan$`F value`[1])
  writeLines(sprintf("Chao1 ANOVA: F=%.3f, P=%.4g", sm_chao$`F value`[1], sm_chao$`Pr(>F)`[1]), "anova_chao1.txt")
  writeLines(sprintf("Shannon ANOVA: F=%.3f, P=%.4g", sm_shan$`F value`[1], sm_shan$`Pr(>F)`[1]), "anova_shannon.txt")
  
  p1 <- ggplot(alpha_df, aes(x=.data[[grp]], y=Chao1, fill=.data[[grp]])) +
    geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.15, alpha=.6, size=1.8) +
    theme_bw(12) + labs(title="Alpha Diversity: Chao1", x=NULL, y="Chao1") +
    annotate("text", x=-Inf, y=Inf, label=lab_chao, hjust=-.1, vjust=1.3, size=4)
  ggsave("alpha_chao1_boxplot.png", p1, width=7, height=5, dpi=200)
  
  p2 <- ggplot(alpha_df, aes(x=.data[[grp]], y=Shannon, fill=.data[[grp]])) +
    geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.15, alpha=.6, size=1.8) +
    theme_bw(12) + labs(title="Alpha Diversity: Shannon", x=NULL, y="Shannon") +
    annotate("text", x=-Inf, y=Inf, label=lab_shan, hjust=-.1, vjust=1.3, size=4)
  ggsave("alpha_shannon_boxplot.png", p2, width=7, height=5, dpi=200)
})

## ===== 02 Beta =====
dir_beta <- mkdir("02_Beta")
with_dir(dir_beta, {
  taxranks <- colnames(tax_table(ps))
  rk <- if ("Genus" %in% taxranks) "Genus" else tail(taxranks,1)
  ps_gen <- suppressWarnings(tax_glom(ps, taxrank=rk, NArm=TRUE))
  ps_rel <- transform_sample_counts(ps_gen, function(x) if (sum(x)==0) x else x/sum(x))
  
  dist_bray <- phyloseq::distance(ps_rel, method="bray")
  ord <- ordinate(ps_rel, method="PCoA", distance=dist_bray)
  
  md <- meta_df %>% tibble::rownames_to_column("Sample") %>% select(Sample, all_of(grp))
  pcoa_df <- as.data.frame(ord$vectors[,1:2]) %>%
    tibble::rownames_to_column("Sample") %>% left_join(md, by="Sample")
  colnames(pcoa_df)[2:3] <- c("PCoA1","PCoA2")
  
  pc1 <- tryCatch(100*ord$values$Relative_eig[1], error=function(e) NA_real_)
  pc2 <- tryCatch(100*ord$values$Relative_eig[2], error=function(e) NA_real_)
  xlab <- if (is.na(pc1)) "PCoA1" else sprintf("PCoA1 (%.1f%%)", pc1)
  ylab <- if (is.na(pc2)) "PCoA2" else sprintf("PCoA2 (%.1f%%)", pc2)
  
  ## adonis2도 formula 필요 → 컬럼명을 안전하게 바꿔 사용
  meta_df_ad <- meta_df
  meta_df_ad$.Group <- meta_df_ad[[grp]]
  ad <- adonis2(dist_bray ~ .Group, data=meta_df_ad, permutations=999)
  Fv <- ad$F[1]; R2 <- ad$R2[1]; Pv <- ad$`Pr(>F)`[1]
  write_tsv(as.data.frame(ad) %>% tibble::rownames_to_column("term"), "permanova_bray_all.tsv")
  
  subt <- sprintf("[PERMANOVA] F-value: %.4f; R-squared: %.5f; p-value: %.3g", Fv, R2, Pv)
  gp <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color=.data[[grp]])) +
    geom_point(size=2, alpha=.9) +
    stat_ellipse(level=.95, linetype=2, linewidth=.5, show.legend=FALSE) +
    theme_bw(12) + labs(title=sprintf("PCoA (Bray, %s)", rk), subtitle=subt, x=xlab, y=ylab)
  ggsave("pcoa_bray_genus.png", gp, width=7, height=5, dpi=200)
  
  ## pairwise PERMANOVA
  groups <- levels(meta_df[[grp]])
  if (length(groups) >= 2) {
    pairs <- t(combn(groups,2)) %>% as.data.frame()
    colnames(pairs) <- c("g1","g2")
    pair_dir <- "pairs"; dir.create(pair_dir, FALSE, TRUE)
    
    res <- purrr::pmap_dfr(pairs, function(g1,g2){
      keep <- rownames(meta_df)[meta_df[[grp]] %in% c(g1,g2)]
      if (length(keep) < 4) return(tibble(g1=g1,g2=g2,F=NA_real_,R2=NA_real_,p=NA_real_,n=length(keep)))
      subD  <- as.dist(as.matrix(dist_bray)[keep, keep])
      subMD <- meta_df[keep, , drop=FALSE]; subMD$.Group <- subMD[[grp]]
      
      ad2 <- adonis2(subD ~ .Group, data=subMD, permutations=999)
      F2 <- ad2$F[1]; R2_ <- ad2$R2[1]; P2 <- ad2$`Pr(>F)`[1]
      
      ps_pair <- prune_samples(keep, ps_rel)
      ord_p <- ordinate(ps_pair, method="PCoA", distance="bray")
      p_df <- as.data.frame(ord_p$vectors[,1:2]) %>%
        tibble::rownames_to_column("Sample") %>%
        left_join(meta_df %>% tibble::rownames_to_column("Sample") %>% select(Sample, all_of(grp)), by="Sample")
      colnames(p_df)[2:3] <- c("PCoA1","PCoA2")
      subt2 <- sprintf("[PERMANOVA] F-value: %.4f; R-squared: %.5f; p-value: %.3g", F2, R2_, P2)
      
      gpp <- ggplot(p_df, aes(PCoA1, PCoA2, color=.data[[grp]])) +
        geom_point(size=2) +
        stat_ellipse(level=.95, linetype=2, linewidth=.5, show.legend=FALSE) +
        theme_bw(12) + labs(title=sprintf("PCoA (Bray, %s vs %s)", g1, g2), subtitle=subt2)
      
      fn <- file.path(pair_dir, sprintf("pcoa_bray_%s_vs_%s.png",
                                        str_replace_all(g1, "[^A-Za-z0-9]+","-"),
                                        str_replace_all(g2, "[^A-Za-z0-9]+","-")))
      ggsave(fn, gpp, width=7, height=5, dpi=200)
      
      tibble(g1=g1,g2=g2,F=F2,R2=R2_,p=P2,n=length(keep))
    })
    res <- res %>% mutate(p_adj = p.adjust(p, method="BH"))
    write_tsv(res, "pairwise_permanova_bray.tsv")
  }
})

## =====================================================================
## 03_Composition (그림 + 모든 랭크 Sam×Grp bar + proportion 테이블 저장)
## =====================================================================
dir_comp <- mkdir("03_Composition")

# ---- A) 기존 Phylum/Genus barplot (샘플×그룹 facet) ----
with_dir(dir_comp, {
  ranks_to_plot <- c("Phylum","Class","Order","Family","Genus","Species")
  for (rk in ranks_to_plot) {
    if (!(rk %in% colnames(tax_table(ps)))) next
    ps_rk  <- suppressWarnings(tax_glom(ps, taxrank = rk, NArm = TRUE))
    ps_rel <- transform_sample_counts(ps_rk, function(x) if (sum(x) == 0) x else x / sum(x))
    
    df <- psmelt(ps_rel) |>
      dplyr::select(Sample, all_of(grp), OTU, Abundance, all_of(rk)) |>
      dplyr::rename(Rank = !!sym(rk))
    
    topN <- if (rk == "Genus") 20 else 10
    top_taxa <- df |>
      dplyr::group_by(Rank) |>
      dplyr::summarise(total = sum(Abundance), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(total)) |>
      dplyr::slice_head(n = topN) |>
      dplyr::pull(Rank)
    
    df <- df |>
      dplyr::mutate(Rank2 = ifelse(Rank %in% top_taxa, Rank, "Others"))
    
    df_sum <- df |>
      dplyr::group_by(Sample, .data[[grp]], Rank2) |>
      dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
    
    p <- ggplot(df_sum, aes(x = Sample, y = Abundance, fill = Rank2)) +
      geom_bar(stat = "identity", width = .95) +
      facet_grid(~ .data[[grp]], scales = "free_x", space = "free_x") +
      theme_bw(12) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
      labs(title = sprintf("Taxonomic composition (%s)", rk), x = "Sample", y = "Relative abundance")
    
    ggsave(sprintf("barplot_%s.png", tolower(rk)), p, width = 14, height = 6, dpi = 200)
  }
})
# 
# # ---- B) MicrobiomeAnalystR: 모든 rank에 대해 PlotTaxaAbundanceBarSamGrp ----
# dir_comp_samgrp <- file.path(dir_comp, "SamGrpBars")
# dir.create(dir_comp_samgrp, showWarnings = FALSE, recursive = TRUE)
# 
# with_dir(dir_comp_samgrp, {
#   ranks_all <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#   
#   run_safely <- function(expr) {
#     tryCatch(force(expr),
#              error = function(e){ message("Skip: ", conditionMessage(e)); invisible(NULL) })
#   }
#   
#   for (rk in ranks_all) {
#     if (!(rk %in% colnames(tax_table(ps)))) {
#       message("Skip SamGrp bar (rank missing): ", rk); next
#     }
#     img_base <- sprintf("taxa_bar_samgrp_%s", tolower(rk))
#     message(">> PlotTaxaAbundanceBarSamGrp: ", rk)
#     run_safely({
#       # (imgName, taxRank, colorFac, facetFac, type, topNum, palette, sumMethod, featNum, legendPos, showOther, fmt)
#       mbSet <<- PlotTaxaAbundanceBarSamGrp(
#         mbSet,
#         img_base,
#         rk,
#         grp,     # color
#         grp,     # facet
#         "barnorm",
#         10,
#         "set3",
#         "sum",
#         10,
#         "bottom",
#         "F",
#         "png"
#       )
#     })
#   }
# })


# ---- B) SamGrpBars를 ggplot으로 직접 생성 (샘플×그룹 & 그룹평균) ----
dir_comp_samgrp <- file.path(dir_comp, "SamGrpBars")
dir.create(dir_comp_samgrp, showWarnings = FALSE, recursive = TRUE)

.clean_taxon <- function(x, rk) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- paste0("Unclassified_", rk)
  x
}

make_samgrp_plots_one_rank <- function(ps, rk, grp, outdir, topN = 10) {
  if (!(rk %in% colnames(tax_table(ps)))) {
    message("Skip SamGrpBars (rank missing): ", rk)
    return(invisible(NULL))
  }
  
  # rank로 집계 후 sample별 상대풍부도
  ps_rk  <- suppressWarnings(tax_glom(ps, taxrank = rk, NArm = TRUE))
  ps_rel <- transform_sample_counts(ps_rk, function(x) if (sum(x) == 0) x else x / sum(x))
  
  df <- psmelt(ps_rel) |>
    dplyr::select(Sample, all_of(grp), OTU, Abundance, all_of(rk)) |>
    dplyr::rename(Group = !!sym(grp), Taxon = !!sym(rk)) |>
    dplyr::mutate(Taxon = .clean_taxon(Taxon, rk))
  
  # 상위 topN taxon만 표기, 나머지 Others
  top_taxa <- df |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(total = sum(Abundance), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::slice_head(n = topN) |>
    dplyr::pull(Taxon)
  
  df <- df |>
    dplyr::mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Others"))
  
  # (1) 샘플별 스택바: x=Sample, facet=Group  → SamGrpBars_sample
  df_samp <- df |>
    dplyr::group_by(Sample, Group, Taxon2) |>
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # 그룹 내 샘플 순서를 Group별로 묶어서 정렬
  sample_order <- df_samp |>
    dplyr::group_by(Group, Sample) |>
    dplyr::summarise(tot = sum(Abundance), .groups = "drop") |>
    dplyr::arrange(Group, dplyr::desc(tot)) |>
    dplyr::pull(Sample)
  
  df_samp$Sample <- factor(df_samp$Sample, levels = unique(sample_order))
  
  p_samp <- ggplot(df_samp, aes(x = Sample, y = Abundance, fill = Taxon2)) +
    geom_col(width = 0.95) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    labs(title = sprintf("Taxonomic composition (Sample × Group) • %s", rk),
         x = "Sample", y = "Relative abundance")
  
  fn1 <- file.path(outdir, sprintf("taxa_bar_SamGrpBars_sample_%s.png", tolower(rk)))
  ggsave(fn1, p_samp, width = 14, height = 6, dpi = 200)
  ggsave(sub("\\.png$", ".pdf", fn1), p_samp, width = 14, height = 6)
  
  # (2) 그룹평균 스택바: x=Group, y=mean(rel)  → SamGrpBars_groupmean
  df_group <- df |>
    dplyr::group_by(Group, Taxon2, Sample) |>
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") |>
    dplyr::group_by(Group, Taxon2) |>
    dplyr::summarise(mean_rel = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  p_grp <- ggplot(df_group, aes(x = Group, y = mean_rel, fill = Taxon2)) +
    geom_col(width = 0.9) +
    theme_bw(12) +
    labs(title = sprintf("Taxonomic composition (Group mean) • %s", rk),
         x = "Group", y = "Mean relative abundance")
  
  fn2 <- file.path(outdir, sprintf("taxa_bar_SamGrpBars_groupmean_%s.png", tolower(rk)))
  ggsave(fn2, p_grp, width = 9, height = 5, dpi = 200)
  ggsave(sub("\\.png$", ".pdf", fn2), p_grp, width = 9, height = 5)
  
  invisible(list(sample_plot = fn1, group_plot = fn2))
}

withr::with_dir(dir_comp_samgrp, {
  ranks_all <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  for (rk in ranks_all) {
    message(">> SamGrpBars ggplot: ", rk)
    try(make_samgrp_plots_one_rank(ps, rk, grp, getwd(), topN = 10), silent = TRUE)
  }
})

# ---- C) 모든 rank에 대해 proportion 테이블 저장 (long / wide / 그룹평균 / 그룹중앙값 / 총합) ----
dir_props <- file.path(dir_comp, "props")
dir.create(dir_props, showWarnings = FALSE, recursive = TRUE)

.clean_taxon <- function(x, rk) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- paste0("Unclassified_", rk)
  x
}

save_rank_props <- function(ps, rk, grp, meta_df, outdir) {
  if (!(rk %in% colnames(tax_table(ps)))) {
    message("Skip props (rank missing): ", rk)
    return(invisible(NULL))
  }
  ps_rk  <- suppressWarnings(tax_glom(ps, taxrank = rk, NArm = TRUE))
  ps_rel <- transform_sample_counts(ps_rk, function(x) if (sum(x) == 0) x else x / sum(x))
  
  # long: Sample, Group, Taxon, Abundance
  df_long <- psmelt(ps_rel) |>
    dplyr::select(Sample, all_of(grp), Abundance, all_of(rk)) |>
    dplyr::rename(Group = !!sym(grp), Taxon = !!sym(rk)) |>
    dplyr::mutate(Taxon = .clean_taxon(Taxon, rk)) |>
    dplyr::arrange(Sample, Taxon)
  
  readr::write_tsv(df_long, file.path(outdir, sprintf("%s_long.tsv", tolower(rk))))
  
  # wide: Sample × Taxon (+Group)
  df_wide <- df_long |>
    dplyr::group_by(Sample, Taxon) |>
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0) |>
    dplyr::left_join(
      meta_df |>
        tibble::rownames_to_column("Sample") |>
        dplyr::select(Sample, all_of(grp)) |>
        dplyr::rename(Group = !!sym(grp)),
      by = "Sample"
    ) |>
    dplyr::relocate(Group, .after = Sample)
  
  readr::write_tsv(df_wide, file.path(outdir, sprintf("%s_wide_samplexTaxon.tsv", tolower(rk))))
  
  # 그룹 평균/중앙값
  df_group_mean <- df_long |>
    dplyr::group_by(Group, Taxon) |>
    dplyr::summarise(mean_prop = mean(Abundance, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Taxon, values_from = mean_prop, values_fill = 0) |>
    dplyr::arrange(Group)
  
  df_group_median <- df_long |>
    dplyr::group_by(Group, Taxon) |>
    dplyr::summarise(median_prop = median(Abundance, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Taxon, values_from = median_prop, values_fill = 0) |>
    dplyr::arrange(Group)
  
  readr::write_tsv(df_group_mean,   file.path(outdir, sprintf("%s_group_mean.tsv",   tolower(rk))))
  readr::write_tsv(df_group_median, file.path(outdir, sprintf("%s_group_median.tsv", tolower(rk))))
  
  # 전체 합 기준 상위 taxa
  df_totals <- df_long |>
    dplyr::group_by(Taxon) |>
    dplyr::summarise(total_prop = sum(Abundance), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total_prop))
  
  readr::write_tsv(df_totals, file.path(outdir, sprintf("%s_taxa_totals.tsv", tolower(rk))))
}

ranks_to_save <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
for (rk in ranks_to_save) {
  message("Saving proportions for rank: ", rk)
  save_rank_props(ps, rk, grp, meta_df, dir_props)
}



cat("== Done ==\nOutput root: ", outdir,
    "\nSubfolders:\n  - 00_QC_Normalization\n  - 01_Alpha\n  - 02_Beta\n  - 03_Composition\n", sep="")
