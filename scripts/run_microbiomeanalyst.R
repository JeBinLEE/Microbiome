#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  req <- c("optparse","data.table","ggplot2","pheatmap","reshape2","vegan","phyloseq")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "))
  lapply(req, function(p) suppressPackageStartupMessages(library(p, character.only=TRUE)))
}))

# ---------- CLI ----------
suppressPackageStartupMessages(library(optparse))
parser <- OptionParser()
parser <- add_option(parser, "--meta", type="character", help="metadata txt (columns: Sample, Group[optional])")
parser <- add_option(parser, "--tax",  type="character", help="taxonomy txt (columns: FeatureID, Taxon)")
parser <- add_option(parser, "--feat", type="character", help="feature-table txt exported from QIIME2")
parser <- add_option(parser, "--out",  type="character", default=".", help="output directory")
opt <- parse_args(parser)
if (is.null(opt$meta) || is.null(opt$tax) || is.null(opt$feat)) {
  stop("Provide --meta, --tax, --feat")
}
dir.create(opt$out, showWarnings=FALSE, recursive=TRUE)

# ---------- Read inputs ----------
suppressPackageStartupMessages(library(data.table))

# metadata
meta <- fread(opt$meta, sep="\t", header=TRUE)
if (!"Sample" %in% names(meta)) stop("metadata must contain 'Sample'")
if (!"Group" %in% names(meta)) meta$Group <- "G1"
meta$Group <- factor(meta$Group)
rownames(meta) <- meta$Sample

# taxonomy
tax <- fread(opt$tax, sep="\t", header=TRUE)
if (!"FeatureID" %in% names(tax)) setnames(tax, names(tax)[1], "FeatureID")
if (!"Taxon" %in% names(tax)) stop("taxonomy needs 'Taxon'")

# feature table (remove commented lines)
strip_hash <- function(fp){
  ln <- readLines(fp, warn=FALSE)
  ln <- ln[!grepl("^\\s*#", ln)]
  tf <- tempfile(fileext=".tsv"); writeLines(ln, tf); tf
}
feat_file <- strip_hash(opt$feat)
ft <- fread(feat_file, sep="\t", header=TRUE, check.names=FALSE)
if ("Feature ID" %in% names(ft)) setnames(ft, "Feature ID", "FeatureID")
if (!"FeatureID" %in% names(ft)) setnames(ft, names(ft)[1], "FeatureID")

ft_ids <- ft$FeatureID
ft_mat <- as.data.frame(ft[, -1, with=FALSE])
rownames(ft_mat) <- ft_ids

# match samples
common_samples <- intersect(colnames(ft_mat), meta$Sample)
if (length(common_samples) < 2) stop("Need >=2 overlapping samples between feature-table and metadata")
ft_mat <- ft_mat[, common_samples, drop=FALSE]
meta2  <- meta[common_samples, , drop=FALSE]

# taxonomy matrix
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
split_tax <- strsplit(tax$Taxon, ";\\s*")
tx <- matrix(NA_character_, nrow=nrow(tax), ncol=length(ranks),
             dimnames=list(tax$FeatureID, ranks))
for (i in seq_len(nrow(tax))) {
  x <- split_tax[[i]]; len <- min(length(x), length(ranks))
  if (len>0) tx[i, 1:len] <- x[1:len]
}
tx <- tx[rownames(ft_mat), , drop=FALSE]

# phyloseq object
suppressPackageStartupMessages(library(phyloseq))
OTU <- otu_table(as.matrix(ft_mat), taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(tx))
SAM <- sample_data(meta2)
ps  <- phyloseq(OTU, TAX, SAM)

# ---------- 0) Library size ----------
libsizes <- colSums(otu_table(ps))
df_lib <- data.frame(Sample=names(libsizes), LibSize=as.numeric(libsizes), Group=meta2$Group)
p_lib <- ggplot2::ggplot(df_lib, ggplot2::aes(x=Sample, y=LibSize, fill=Group)) +
  ggplot2::geom_col() + ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, hjust=1))
ggplot2::ggsave(file.path(opt$out,"norm_libsizes_0.png"), p_lib, width=10, height=4, dpi=150)

# ---------- 1) Alpha (Shannon) ----------
X <- t(as.matrix(otu_table(ps))) # samples x taxa
sh <- vegan::diversity(X, index="shannon")
alpha_tab <- data.frame(Sample=names(sh), shannon=as.numeric(sh), Group=meta2$Group)
data.table::fwrite(alpha_tab, file.path(opt$out,"alpha_diversity.tsv"), sep="\t")
p_a <- ggplot2::ggplot(alpha_tab, ggplot2::aes(x=Group, y=shannon, fill=Group)) +
  ggplot2::geom_boxplot(width=0.6) + ggplot2::geom_jitter(width=0.1, alpha=0.6) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(opt$out,"alpha_shannon_boxplot.png"), p_a, width=5, height=4, dpi=150)

# ---------- 2) Beta (Bray) + PCoA + PERMANOVA ----------
X_rel <- sweep(X, 1, pmax(rowSums(X),1), "/")
d_bray <- vegan::vegdist(X_rel, method="bray")
coord  <- cmdscale(d_bray, k=2, eig=TRUE)
df_pcoa <- data.frame(Sample=rownames(coord$points),
                      PC1=coord$points[,1], PC2=coord$points[,2],
                      Group=meta2$Group)
p_p <- ggplot2::ggplot(df_pcoa, ggplot2::aes(PC1, PC2, color=Group)) +
  ggplot2::geom_point(size=3) + ggplot2::theme_bw()
ggplot2::ggsave(file.path(opt$out,"beta_bray_pcoa.png"), p_p, width=5, height=4, dpi=150)
perm <- vegan::adonis2(d_bray ~ Group, data=meta2)
sink(file.path(opt$out,"permanova_bray.tsv")); print(perm); sink()

# ---------- 3) "LEfSe 대체": Kruskal–Wallis + epsilon^2 ----------
target_rank <- if ("Genus" %in% colnames(tax_table(ps))) "Genus"
else if ("Family" %in% colnames(tax_table(ps))) "Family" else NA

if (!is.na(target_rank)) {
  ps_rank <- tax_glom(ps, taxrank=target_rank, NArm=TRUE)
  ps_rel  <- transform_sample_counts(ps_rank, function(x) if (sum(x)==0) x else x/sum(x))
} else {
  ps_rel <- transform_sample_counts(ps, function(x) if (sum(x)==0) x else x/sum(x))
  target_rank <- "ASV"
}

sdf <- as(sample_data(ps_rel), "data.frame")
if (nlevels(sdf$Group) >= 2) {
  abund <- as.data.frame(t(otu_table(ps_rel))) # samples x taxa
  dat <- cbind(Group=sdf$Group, abund)
  k <- nlevels(dat$Group); n <- nrow(dat)
  res <- lapply(colnames(dat)[-1], function(t){
    x <- dat[[t]]
    kw <- try(suppressWarnings(kruskal.test(x ~ dat$Group)), silent=TRUE)
    if (inherits(kw, "try-error")) return(NULL)
    H <- as.numeric(kw$statistic); p <- as.numeric(kw$p.value)
    eps2 <- if (n > k) (H - k + 1) / (n - k) else NA_real_
    data.frame(taxon=t, kw_stat=H, pval=p, eps2=eps2, stringsAsFactors=FALSE)
  })
  res <- data.table::rbindlist(res, fill=TRUE)
  res$padj <- p.adjust(res$pval, method="fdr")
  res <- res[order(res$padj, -res$eps2), ]
  data.table::fwrite(res, file.path(opt$out,"lefse_markers.tsv"), sep="\t")
  
  topn <- head(res[!is.na(res$eps2), ], 20)
  out_plot <- file.path(opt$out, "lefse_LDA_barplot.png")
  if (nrow(topn) > 0) {
    topn$taxon <- factor(topn$taxon, levels=rev(topn$taxon))
    gp <- ggplot2::ggplot(topn, ggplot2::aes(x=taxon, y=eps2)) +
      ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw() +
      ggplot2::labs(title=sprintf("Markers by effect size (%s, epsilon^2)", target_rank),
                    x=target_rank, y="Effect size (epsilon^2)")
    ggplot2::ggsave(out_plot, gp, width=8, height=6, dpi=150)
  } else {
    file.create(out_plot)
  }
} else {
  data.table::fwrite(data.table::data.table(message="Only one group; skip KW/effect size"),
                     file.path(opt$out,"lefse_markers.tsv"), sep="\t")
  file.create(file.path(opt$out,"lefse_LDA_barplot.png"))
}

# ---------- 4) Correlation (Spearman): numeric metadata vs top 30 taxa ----------
sdf <- as(sample_data(ps_rel), "data.frame")
num_cols <- names(Filter(is.numeric, sdf))
out_corr <- file.path(opt$out,"correlations_spearman.tsv")
out_heat <- file.path(opt$out,"corr_heatmap_spearman.png")

if (length(num_cols) == 0) {
  data.table::fwrite(data.table::data.table(message="No numeric metadata; skip correlations"),
                     out_corr, sep="\t")
  file.create(out_heat)
} else {
  A <- as.data.frame(t(otu_table(ps_rel)))
  means <- sort(colMeans(A), decreasing=TRUE)
  top_taxa <- names(head(means, 30))
  A <- A[, top_taxa, drop=FALSE]
  cors <- list()
  for (m in num_cols) {
    for (t in colnames(A)) {
      ct <- suppressWarnings(cor.test(sdf[[m]], A[[t]], method="spearman"))
      cors[[length(cors)+1]] <- data.frame(meta=m, taxon=t, rho=ct$estimate, pval=ct$p.value)
    }
  }
  C <- data.table::rbindlist(cors)
  C$padj <- p.adjust(C$pval, "fdr")
  data.table::fwrite(C, out_corr, sep="\t")
  
  M <- reshape2::dcast(C, meta ~ taxon, value.var="rho", fun.aggregate=mean)
  rn <- M$meta; M$meta <- NULL; mat <- as.matrix(M); rownames(mat) <- rn
  if (nrow(mat) >= 1 && ncol(mat) >= 1) {
    pheatmap::pheatmap(mat, filename=out_heat, width=10, height=6)
  } else {
    file.create(out_heat)
  }
}

cat("Done.\n")
