#!/usr/bin/env Rscript
.libPaths("/home/ljb/R/x86_64-pc-linux-gnu-library/4.4")

# Microbiome analysis pipeline (R) with sweep-able parameters.
# - Inputs: feature table (rows=feature, cols=samples), taxonomy table (rows=feature), metadata table (rows=samples)
# - Outputs: standardized folders under --outdir:
#     00_QC_Normalization/
#     01_Alpha/
#     02_Beta/
#     03_Composition/
#
# Notes:
# - Designed to be driven by Snakefile3 sweeps.
# - Does NOT depend on the MicrobiomeAnalyst web app; it reproduces common steps using phyloseq/vegan/ggplot2.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(phyloseq)
  library(vegan)
})

# -------------------------------------------------
# Fail-safe helpers: keep snakemake running, but leave a clear error marker
# -------------------------------------------------
write_error_and_quit <- function(outdir, msg) {
  dir_create(outdir)
  writeLines(as.character(msg), file.path(outdir, "PIPELINE_ERROR.txt"))
  message(msg)
  quit(save = "no", status = 0)
}

# -----------------------
# CLI
# -----------------------
opt_list <- list(
  make_option("--feat_fp", type = "character", help = "Feature table (MA format)"),
  make_option("--tax_fp", type = "character", help = "Taxonomy table (MA format)"),
  make_option("--meta_fp", type = "character", help = "Metadata table (MA format)"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--group", type = "character", default = "Group1", help = "Grouping column in metadata"),

  # filtering
  make_option("--filter_enabled",
    type = "character", default = "true",
    help = "true/false: apply abundance + prevalence + IQR filters"
  ),
  make_option("--min_total_count", type = "double", default = 4, help = "Min total count per feature"),
  make_option("--prevalence", type = "double", default = 0.10, help = "Min prevalence fraction (0-1)"),
  make_option("--iqr_remove", type = "double", default = 0.10, help = "Remove lowest IQR fraction (0-1)"),

  # normalization
  make_option("--scale",
    type = "character", default = "rarefy",
    help = "none|tss|rarefy"
  ),
  make_option("--rarefy_depth",
    type = "integer", default = 0,
    help = "Depth for rarefaction (0 => use min library size)"
  ),
  make_option("--rarefy_q",
    type = "double", default = 0.0,
    help = "Quantile for rarefaction depth if rarefy_depth=0 (0-1). Ignored if rarefy_depth>0"
  ),
  make_option("--transform",
    type = "character", default = "none",
    help = "none|log|clr"
  ),
  make_option("--pseudocount", type = "double", default = 1, help = "Pseudocount for log/clr transforms"),
  make_option("--seed", type = "integer", default = 22, help = "Random seed"),

  # alpha/beta selections (comma-separated)
  make_option("--alpha_ranks",
    type = "character", default = "ASV,Genus,Family",
    help = "Comma-separated ranks for alpha diversity"
  ),
  make_option("--alpha_measures",
    type = "character", default = "Observed,Chao1,Shannon,Simpson",
    help = "Comma-separated alpha measures"
  ),
  make_option("--alpha_stats",
    type = "character", default = "kruskal,wilcoxon,anova",
    help = "Comma-separated alpha stats methods"
  ),
  make_option("--beta_ranks",
    type = "character", default = "Genus,Family",
    help = "Comma-separated ranks for beta diversity"
  ),
  make_option("--beta_distances",
    type = "character", default = "bray,jaccard,euclidean",
    help = "Comma-separated beta distances"
  ),
  make_option("--beta_ordinations",
    type = "character", default = "PCoA,NMDS",
    help = "Comma-separated ordinations"
  ),
  make_option("--beta_stats",
    type = "character", default = "permanova,anosim",
    help = "Comma-separated beta stats methods"
  ),
  make_option("--permutations", type = "integer", default = 999, help = "Permutations for PERMANOVA/ANOSIM")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# -----------------------
# helpers
# -----------------------
dir_create <- function(p) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || x == "") {
    return(character())
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

is_true <- function(x) tolower(trimws(x)) %in% c("1", "true", "t", "yes", "y")

as_frac <- function(x) {
  if (is.null(x) || is.na(x)) {
    return(x)
  }
  if (is.numeric(x) && length(x) == 1 && x > 1) {
    return(x / 100)
  }
  return(x)
}

safe_tax_glom <- function(ps, rk) {
  if (rk == "ASV") {
    return(ps)
  }
  tt <- tryCatch(tax_table(ps), error = function(e) NULL)
  if (is.null(tt)) {
    return(NULL)
  }
  cn <- colnames(tt)
  if (is.null(cn) || !(rk %in% cn)) {
    return(NULL)
  }
  out <- tryCatch(tax_glom(ps, taxrank = rk, NArm = FALSE), error = function(e) NULL)
  return(out)
}

tidy_kruskal <- function(test) {
  data.frame(statistic = unname(test$statistic), df = unname(test$parameter), p.value = test$p.value, check.names = FALSE)
}
tidy_wilcox <- function(test) {
  data.frame(statistic = unname(test$statistic), p.value = test$p.value, check.names = FALSE)
}
tidy_anova <- function(fit) {
  sm <- summary(fit)[[1]]
  sm <- data.frame(term = rownames(sm), sm, check.names = FALSE)
  rownames(sm) <- NULL
  sm
}


# For taxonomy ranks (case-insensitive match)
normalize_rank <- function(r) {
  r2 <- trimws(r)
  if (tolower(r2) == "asv") {
    return("ASV")
  }
  # Common taxonomy column names in MA taxonomy_for_MA.txt:
  # Kingdom Phylum Class Order Family Genus Species
  # Keep as-is
  r2
}

# -----------------------
# I/O
# -----------------------
set.seed(opt$seed)

outdir <- opt$outdir
dir_create(outdir)

qc_dir <- file.path(outdir, "00_QC_Normalization")
alpha_dir <- file.path(outdir, "01_Alpha")
beta_dir <- file.path(outdir, "02_Beta")
comp_dir <- file.path(outdir, "03_Composition")
dir_create(qc_dir)
dir_create(alpha_dir)
dir_create(beta_dir)
dir_create(comp_dir)

# Read tables
feat <- fread(opt$feat_fp, data.table = FALSE)
stopifnot(ncol(feat) >= 3)
rownames(feat) <- feat[[1]]
feat[[1]] <- NULL
feat_mat <- as.matrix(feat)
mode(feat_mat) <- "numeric"

tax <- fread(opt$tax_fp, data.table = FALSE)
rownames(tax) <- tax[[1]]
tax[[1]] <- NULL

meta <- fread(opt$meta_fp, data.table = FALSE)
if (!("#NAME" %in% colnames(meta))) {
  # allow "Sample" only
  if ("Sample" %in% colnames(meta)) meta$`#NAME` <- meta$Sample
}
stopifnot("#NAME" %in% colnames(meta))
rownames(meta) <- meta$`#NAME`

# align samples
common_samples <- intersect(colnames(feat_mat), rownames(meta))
if (length(common_samples) < 2) stop("Too few matched samples between feature table and metadata.")
feat_mat <- feat_mat[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# group column
if (!(opt$group %in% colnames(meta))) {
  stop(sprintf(
    "Group column '%s' not found in metadata columns: %s",
    opt$group, paste(colnames(meta), collapse = ", ")
  ))
}
meta[[opt$group]] <- as.factor(meta[[opt$group]])

# align taxonomy rows
common_taxa <- intersect(rownames(feat_mat), rownames(tax))
feat_mat <- feat_mat[common_taxa, , drop = FALSE]
tax <- tax[common_taxa, , drop = FALSE]

# Build phyloseq object
OTU <- otu_table(feat_mat, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
SAM <- sample_data(as.data.frame(meta))
ps0 <- phyloseq(OTU, TAX, SAM)

# -----------------------
# Filtering
# -----------------------
n_taxa0 <- ntaxa(ps0)
n_samp0 <- nsamples(ps0)

ps_f <- ps0
if (is_true(opt$filter_enabled)) {
  otu <- as(otu_table(ps_f), "matrix")
  # abundance + prevalence (MicrobiomeAnalyst-style)
  # - keep features that have count >= count_cut in at least `prevalence` fraction of samples
  # - also require total count >= count_cut (light safeguard)
  count_cut <- opt$min_total_count
  total <- rowSums(otu)
  prev <- rowSums(otu >= count_cut) / ncol(otu)

  keep1 <- (total >= count_cut) & (prev >= opt$prevalence)
  n_taxa_after_lowcount <- sum(keep1, na.rm = TRUE)
  ps_f <- prune_taxa(keep1, ps_f)

  # IQR filter (remove lowest fraction)
  otu2 <- as(otu_table(ps_f), "matrix")
  if (nrow(otu2) > 0 && opt$iqr_remove > 0) {
    iqr <- apply(otu2, 1, IQR)
    n_rm <- floor(length(iqr) * opt$iqr_remove)
    if (n_rm > 0) {
      rm_ids <- names(sort(iqr, decreasing = FALSE))[seq_len(n_rm)]
      ps_f <- prune_taxa(!(taxa_names(ps_f) %in% rm_ids), ps_f)
    }
  }
}

# Drop samples with zero total counts after filtering (common cause of downstream rarefaction failures)
libs_tmp <- tryCatch(phyloseq::sample_sums(ps_f), error = function(e) NULL)
if (!is.null(libs_tmp)) {
  zero_samps <- names(libs_tmp)[libs_tmp <= 0]
  if (length(zero_samps)) {
    message("[filter] Dropping ", length(zero_samps), " samples with zero total counts after filtering")
    ps_f <- prune_samples(setdiff(sample_names(ps_f), zero_samps), ps_f)
  }
}

# Hard stop conditions (write error marker and exit 0 so snakemake can continue)
if (nsamples(ps_f) < 2) {
  write_error_and_quit(opt$outdir, paste0("Too few samples remain after filtering: ", nsamples(ps_f)))
}
if (ntaxa(ps_f) < 1) {
  write_error_and_quit(opt$outdir, "No features remain after filtering")
}


# Filter breakdown summary
n_taxa_after_iqr <- ntaxa(ps_f)
filter_details <- data.frame(
  stage = c("before", "after_lowcount_prevalence", "after_iqr"),
  n_samples = c(n_samp0, n_samp0, nsamples(ps_f)),
  n_features = c(n_taxa0, if (exists("n_taxa_after_lowcount")) n_taxa_after_lowcount else NA, n_taxa_after_iqr),
  count_cut = opt$min_total_count,
  prevalence = opt$prevalence,
  iqr_remove = opt$iqr_remove,
  stringsAsFactors = FALSE
)
fwrite(filter_details, file.path(qc_dir, "filter_details.tsv"), sep = "\t")

n_taxaF <- ntaxa(ps_f)
n_sampF <- nsamples(ps_f)

# Save filtered tables for downstream
feat_f <- as(otu_table(ps_f), "matrix")
tax_f <- as(tax_table(ps_f), "matrix")
meta_f <- as(sample_data(ps_f), "data.frame")

fwrite(
  data.frame("#NAME" = rownames(feat_f), feat_f, check.names = FALSE),
  file.path(qc_dir, "feature_table_filtered.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  data.frame("#TAXONOMY" = rownames(tax_f), tax_f, check.names = FALSE),
  file.path(qc_dir, "taxonomy_filtered.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  data.frame("#NAME" = rownames(meta_f), meta_f, check.names = FALSE),
  file.path(qc_dir, "metadata_filtered.tsv"),
  sep = "\t", quote = FALSE
)

qc_sum <- data.frame(
  metric = c("n_samples", "n_features"),
  before = c(n_samp0, n_taxa0),
  after  = c(n_sampF, n_taxaF)
)
fwrite(qc_sum, file.path(qc_dir, "qc_summary.tsv"), sep = "\t")

qc_long <- data.frame(
  key = c(
    "n_samples_before", "n_samples_after",
    "n_features_before", "n_features_after_lowcount_prevalence", "n_features_after_iqr"
  ),
  value = c(
    n_samp0, n_sampF,
    n_taxa0,
    if (exists("n_taxa_after_lowcount")) n_taxa_after_lowcount else NA,
    if (exists("n_taxa_after_iqr")) n_taxa_after_iqr else n_taxaF
  )
)
fwrite(qc_long, file.path(qc_dir, "qc_summary_long.tsv"), sep = "\t")

# Library size plot (after filtering)
lib <- colSums(feat_f)
lib_df <- data.frame(sample = names(lib), libsize = as.numeric(lib), group = meta_f[[opt$group]])
p_lib <- ggplot(lib_df, aes(x = group, y = libsize)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0) +
  labs(title = "Library size after filtering", x = opt$group, y = "Reads") +
  theme_bw()
ggsave(file.path(qc_dir, "libsize_boxplot.pdf"), p_lib, width = 7, height = 4)
fwrite(lib_df, file.path(qc_dir, "libsize.tsv"), sep = "\t")

# -----------------------
# Normalization / transform
# -----------------------
ps_scaled <- ps_f
scale <- tolower(opt$scale)
if (scale == "rarefy") {
  depth <- opt$rarefy_depth
  if (depth <= 0) {
    if (!is.na(opt$rarefy_q) && opt$rarefy_q > 0 && opt$rarefy_q <= 1) {
      depth <- floor(as.numeric(stats::quantile(lib, probs = opt$rarefy_q, na.rm = TRUE)))
    } else {
      depth <- min(lib)
    }
  }
  depth <- min(depth, min(lib, na.rm = TRUE))
  ps_scaled <- rarefy_even_depth(ps_scaled, sample.size = depth, rngseed = opt$seed, replace = FALSE, verbose = FALSE)
} else if (scale == "tss") {
  otu <- as(otu_table(ps_scaled), "matrix")
  otu <- sweep(otu, 2, colSums(otu), "/")
  otu[is.na(otu)] <- 0
  otu_table(ps_scaled) <- otu_table(otu, taxa_are_rows = TRUE)
} else if (scale == "none") {
  # do nothing
} else {
  stop("Unsupported --scale. Use none|tss|rarefy")
}

# Store scaled (positive) table for Bray/Jaccard
mat_pos <- as(otu_table(ps_scaled), "matrix")

# transform
ps_final <- ps_scaled
tr <- tolower(opt$transform)
if (tr == "log") {
  otu <- as(otu_table(ps_final), "matrix")
  otu <- log(otu + opt$pseudocount)
  otu_table(ps_final) <- otu_table(otu, taxa_are_rows = TRUE)
} else if (tr == "clr") {
  otu <- as(otu_table(ps_final), "matrix")
  otu <- otu + opt$pseudocount
  logx <- log(otu)
  clr <- sweep(logx, 2, colMeans(logx), "-")
  otu_table(ps_final) <- otu_table(clr, taxa_are_rows = TRUE)
} else if (tr == "none") {
  # do nothing
} else {
  stop("Unsupported --transform. Use none|log|clr")
}

# Save normalized table (final)
feat_norm <- as(otu_table(ps_final), "matrix")
fwrite(
  data.frame("#NAME" = rownames(feat_norm), feat_norm, check.names = FALSE),
  file.path(qc_dir, "feature_table_normalized.tsv"),
  sep = "\t", quote = FALSE
)

# -----------------------
# Alpha diversity
# -----------------------
alpha_ranks <- vapply(split_csv(opt$alpha_ranks), normalize_rank, character(1))
alpha_meas <- split_csv(opt$alpha_measures)
alpha_stats <- split_csv(opt$alpha_stats)

for (rk in alpha_ranks) {
  rk_dir <- file.path(alpha_dir, rk)
  dir_create(rk_dir)

  ps_rk <- ps_scaled
  if (rk != "ASV") {
    ps_rk2 <- safe_tax_glom(ps_rk, rk)
    if (is.null(ps_rk2)) {
      message("Skipping alpha rank (not present in taxonomy): ", rk)
      next
    }
    ps_rk <- ps_rk2
  }

  # estimate_richness uses counts -> use ps_scaled
  # TSS/CLR produce non-integer data, which will fail estimate_richness
  # Wrap in tryCatch to skip gracefully
  rich <- tryCatch(
    {
      estimate_richness(ps_rk, measures = alpha_meas)
    },
    error = function(e) {
      message(
        "Warning: Alpha diversity calculation failed for rank ", rk,
        " (likely due to non-integer data from TSS/CLR normalization). Skipping."
      )
      message("Error message: ", e$message)
      return(NULL)
    }
  )

  if (is.null(rich)) {
    # Write a note file explaining why alpha was skipped
    writeLines(
      paste0(
        "Alpha diversity skipped for this normalization: ",
        "estimate_richness requires integer count data, but TSS/CLR produces proportions."
      ),
      file.path(rk_dir, "ALPHA_SKIPPED_NOTE.txt")
    )
    next
  }

  rich$sample <- rownames(rich)
  md <- as(sample_data(ps_rk), "data.frame")
  md$sample <- rownames(md)
  df <- left_join(rich, md, by = "sample")

  fwrite(df, file.path(rk_dir, "alpha_diversity.tsv"), sep = "\t")

  # plots + stats per measure
  for (m in alpha_meas) {
    if (!(m %in% colnames(df))) next
    p <- ggplot(df, aes(x = .data[[opt$group]], y = .data[[m]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, height = 0) +
      labs(title = paste0("Alpha: ", rk, " | ", m), x = opt$group, y = m) +
      theme_bw()
    ggsave(file.path(rk_dir, paste0("alpha_", m, ".pdf")), p, width = 7, height = 4)

    res <- list()
    g0 <- df[[opt$group]]
    y0 <- df[[m]]

    ok <- !is.na(g0) & !is.na(y0)
    g <- factor(g0[ok])
    y <- y0[ok]

    if ("kruskal" %in% alpha_stats) {
      if (nlevels(g) >= 2) {
        res$kruskal <- tidy_kruskal(kruskal.test(y ~ g))
      }
    }
    if ("anova" %in% alpha_stats) {
      if (nlevels(g) >= 2) {
        res$anova <- tidy_anova(aov(y ~ g))
      }
    }
    if ("wilcoxon" %in% alpha_stats) {
      if (nlevels(g) == 2) {
        res$wilcoxon <- tidy_wilcox(wilcox.test(y ~ g))
      }
    }

    # write stats
    if (length(res) > 0) {
      for (k in names(res)) {
        fwrite(res[[k]], file.path(rk_dir, paste0("alpha_", m, "_", k, ".tsv")), sep = "\t")
      }
    }
  }
}

# -----------------------
# Beta diversity
# -----------------------
beta_ranks <- vapply(split_csv(opt$beta_ranks), normalize_rank, character(1))
beta_dist <- tolower(split_csv(opt$beta_distances))
beta_ord <- toupper(split_csv(opt$beta_ordinations))
beta_stats <- tolower(split_csv(opt$beta_stats))

# helper to get matrix per rank for positive (pos) and final (final)
get_rank_mats <- function(rk) {
  ps_pos <- ps_scaled
  ps_fin <- ps_final
  if (rk != "ASV") {
    ps_pos2 <- safe_tax_glom(ps_pos, rk)
    ps_fin2 <- safe_tax_glom(ps_fin, rk)
    if (is.null(ps_pos2) || is.null(ps_fin2)) {
      return(NULL)
    }
    ps_pos <- ps_pos2
    ps_fin <- ps_fin2
  }
  list(
    pos = t(as(otu_table(ps_pos), "matrix")), # samples x features
    fin = t(as(otu_table(ps_fin), "matrix")), # samples x features
    meta = as(sample_data(ps_pos), "data.frame")
  )
}

for (rk in beta_ranks) {
  rk_dir <- file.path(beta_dir, rk)
  dir_create(rk_dir)

  mats <- get_rank_mats(rk)
  if (is.null(mats)) {
    message("Skipping beta rank (not present in taxonomy): ", rk)
    next
  }
  md <- mats$meta
  grp <- factor(md[[opt$group]])

  for (d in beta_dist) {
    x <- if (d %in% c("bray", "jaccard")) mats$pos else mats$fin
    dist <- vegan::vegdist(x, method = d)

    # stats
    if ("permanova" %in% beta_stats) {
      ad <- vegan::adonis2(dist ~ grp, permutations = opt$permutations)
      fwrite(as.data.frame(ad), file.path(rk_dir, paste0("beta_", d, "_permanova.tsv")), sep = "\t")
    }
    if ("anosim" %in% beta_stats) {
      # ANOSIM
      an <- tryCatch(
        {
          vegan::anosim(dist, grp, permutations = opt$permutations)
        },
        error = function(e) {
          message("Skipping ANOSIM due to error: ", e$message)
          return(NULL)
        }
      )
      if (!is.null(an)) {
        fwrite(data.frame(statistic = an$statistic, signif = an$signif, permutations = an$permutations),
          file.path(rk_dir, paste0("beta_", d, "_anosim.tsv")),
          sep = "\t"
        )
      }
    }

    # ordinations
    for (ord in beta_ord) {
      ord_dir <- file.path(rk_dir, paste0(d, "_", ord))
      dir_create(ord_dir)

      if (ord == "PCOA") {
        co <- cmdscale(dist, k = 2, eig = TRUE)
        df <- data.frame(sample = rownames(co$points), PC1 = co$points[, 1], PC2 = co$points[, 2])
      } else if (ord == "NMDS") {
        nmds <- vegan::metaMDS(dist, k = 2, trymax = 50, autotransform = FALSE, trace = FALSE)
        df <- data.frame(sample = rownames(nmds$points), NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2])
      } else {
        next
      }

      df <- left_join(df, md %>% mutate(sample = rownames(md)), by = "sample")
      fwrite(df, file.path(ord_dir, "ordination_points.tsv"), sep = "\t")

      x1 <- colnames(df)[2]
      x2 <- colnames(df)[3]
      p <- ggplot(df, aes(x = .data[[x1]], y = .data[[x2]], color = .data[[opt$group]])) +
        geom_point(size = 2, alpha = 0.9) +
        labs(title = paste0("Beta: ", rk, " | ", d, " | ", ord), color = opt$group) +
        theme_bw()
      ggsave(file.path(ord_dir, "ordination.pdf"), p, width = 6, height = 5)
    }
  }
}

# -----------------------
# Composition (stacked bar, top taxa)
# -----------------------
comp_ranks <- unique(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
for (rk in comp_ranks) {
  rk_dir <- file.path(comp_dir, rk)
  dir_create(rk_dir)
  ps_rk <- safe_tax_glom(ps_scaled, rk)
  if (is.null(ps_rk)) {
    message("Skipping composition rank (not present in taxonomy): ", rk)
    next
  }
  otu <- as(otu_table(ps_rk), "matrix")
  otu_rel <- sweep(otu, 2, colSums(otu), "/")
  otu_rel[is.na(otu_rel)] <- 0

  taxm <- as(tax_table(ps_rk), "matrix")
  taxa_label <- as.character(taxm[, rk])
  taxa_label[is.na(taxa_label) | taxa_label == ""] <- "Unassigned"
  rownames(otu_rel) <- make.unique(taxa_label)

  # top 20 taxa overall
  topN <- 20
  totals <- rowSums(otu_rel)
  top_taxa <- names(sort(totals, decreasing = TRUE))[seq_len(min(topN, length(totals)))]
  other <- setdiff(rownames(otu_rel), top_taxa)

  otu_rel2 <- otu_rel
  if (length(other) > 0) {
    other_row <- colSums(otu_rel2[other, , drop = FALSE])
    otu_rel2 <- otu_rel2[top_taxa, , drop = FALSE]
    otu_rel2 <- rbind(otu_rel2, Other = other_row)
  }

  # Long format dataframe for sample-level plots
  df <- as.data.frame(otu_rel2) %>%
    {
      data.frame(Taxon = rownames(.), ., check.names = FALSE)
    } %>%
    pivot_longer(-Taxon, names_to = "sample", values_to = "abundance") %>%
    left_join(as(sample_data(ps_scaled), "data.frame") %>% mutate(sample = rownames(.)), by = "sample")

  fwrite(df, file.path(rk_dir, "composition_long.tsv"), sep = "\t")

  # 1) Sample-level stacked bar (original, renamed)
  p_sample <- ggplot(df, aes(x = sample, y = abundance, fill = Taxon)) +
    geom_col(width = 0.9) +
    facet_grid(~ .data[[opt$group]], scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.spacing = unit(0.2, "lines")) +
    labs(title = paste0("Composition by Sample: ", rk), x = "Sample", y = "Relative abundance")
  ggsave(file.path(rk_dir, "stacked_bar_by_sample.pdf"), p_sample, width = 12, height = 4)

  # 2) Group-averaged stacked bar
  df_group <- df %>%
    group_by(.data[[opt$group]], Taxon) %>%
    summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop")

  # Order taxa by total abundance for consistent stacking
  taxa_order <- df_group %>%
    group_by(Taxon) %>%
    summarise(total = sum(abundance)) %>%
    arrange(desc(total)) %>%
    pull(Taxon)
  df_group$Taxon <- factor(df_group$Taxon, levels = rev(taxa_order))

  p_group <- ggplot(df_group, aes(x = .data[[opt$group]], y = abundance, fill = Taxon)) +
    geom_col(width = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = paste0("Mean Composition by Group: ", rk), x = opt$group, y = "Relative abundance")
  ggsave(file.path(rk_dir, "stacked_bar_by_group.pdf"), p_group, width = 8, height = 5)

  # Save group-level summary table
  df_group_wide <- df_group %>%
    pivot_wider(names_from = Taxon, values_from = abundance, values_fill = 0)
  fwrite(df_group_wide, file.path(rk_dir, "composition_by_group.tsv"), sep = "\t")

  # 3) Area plot by group (for visualizing composition trends)
  # Reorder for area plot
  df_group_area <- df_group
  df_group_area$Taxon <- factor(df_group_area$Taxon, levels = taxa_order)

  p_area <- ggplot(df_group_area, aes(x = .data[[opt$group]], y = abundance, fill = Taxon, group = Taxon)) +
    geom_area(position = "stack", alpha = 0.85) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = paste0("Composition Area Plot: ", rk), x = opt$group, y = "Relative abundance")
  ggsave(file.path(rk_dir, "area_plot_by_group.pdf"), p_area, width = 8, height = 5)
}

message("DONE: ", outdir)
