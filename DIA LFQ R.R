library("tidyverse")
library("readxl")
library("limma")
library("ggrepel")
library("ggsci")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("writexl")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

need <- c("enrichplot", "DOSE", "GOSemSim", "clusterProfiler")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) BiocManager::install(miss)

library("clusterProfiler")
library("GOSemSim")
library("enrichplot")


# Output setup - folders for exported figures ---------------------------------

fig_root <- "figures"
fig_dirs <- c(qc          = file.path(fig_root, "01_QC"),
              pca         = file.path(fig_root, "02_PCA"),
              volcano     = file.path(fig_root, "03_Volcano"),
              go          = file.path(fig_root, "04_GO"),
              enrichment  = file.path(fig_root, "05_Enrichment"),
              tables      = file.path(fig_root, "06_Tables"),
              concordance = file.path(fig_root, "07_Concordance"))

invisible(lapply(fig_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

## saving figure to a folder as PDF
save_fig <- function(p, folder, name, w = 7, h = 5) {
  if (is.null(p)) { message("Nothing to save for: ", name); return(invisible(NULL)) }
  ggsave(file.path(fig_dirs[[folder]], paste0(name, ".pdf")),
         plot = p, width = w, height = h, device = cairo_pdf)
  print(p)
  invisible(p)
}


# 1. Load the data, make sure your directory is correct. ----------------------

dat <- read_excel("DIA data/260901_AL_ASK_ASK1_proteins.xlsx",
                  sheet = "260901_AL_ASK_ASK1.pg_matrix") %>% as.data.frame()

df <- dat

qcols <- grep("^[A-D]\\d+_", names(df), value = TRUE)

length(qcols) # just to confirm the columns are all counted for.....


# 2. Peptide-count filter -----------------------------------------------------

# NOTE: the >= 2 sequence filter is DISABLED; single-peptide groups are kept.
# Only rows with a missing Protein.Group are dropped (needed as rownames).
df <- df[!is.na(df$Protein.Group), ]
nrow(df)

# ---- To re-enable the two-peptide filter, comment out the line above and
# ---- uncomment the block below:
#
# df <- df[!is.na(df$Protein.Group) &
#            !is.na(df$N.Sequences) &
#            as.numeric(df$N.Sequences) >= 2, ]
# nrow(df)

# -----------------------------------------------------------------------------
# 2b. Remove contaminant protein groups (Cont_ prefix from the search database)
# -----------------------------------------------------------------------------
n_before <- nrow(df)
is_contaminant <- vapply(strsplit(df$Protein.Group, ";"), 
                         function(a) all(grepl("^(Cont_|cRAP-)", a)), logical(1))
df <- df[!is_contaminant, ]



# 3. Build the intensity matrix, log2 transform -------------------------------

# Explicit coercion: data.matrix() would silently turn a character column into
# factor codes instead of intensities, with no error.
num <- lapply(df[qcols], function(x) suppressWarnings(as.numeric(as.character(x))))
bad <- vapply(seq_along(num),
              function(i) sum(is.na(num[[i]]) & !is.na(df[[qcols[i]]])), integer(1))
if (any(bad > 0))
  warning("Non-numeric entries coerced to NA in: ",
          paste(sprintf("%s (%d)", qcols[bad > 0], bad[bad > 0]), collapse = ", "))

mat <- as.matrix(as.data.frame(num))
colnames(mat) <- qcols
mat[mat == 0] <- NA

# Rownames must be unique: match() returns only the first hit, so duplicates
# would silently inherit the wrong annotation.
if (any(duplicated(df$Protein.Group)))
  stop("Duplicated Protein.Group entries: ",
       paste(unique(df$Protein.Group[duplicated(df$Protein.Group)]),
             collapse = ", "))

rownames(mat) <- df$Protein.Group
storage.mode(mat) <- "double"


log_mat <- log2(mat)

# -----------------------------------------------------------------------------
#Missing-value filter: drop proteins quantified in too few samples.
#Keep a protein with >= min_valid non-missing values in >= 1 group
# -----------------------------------------------------------------------------


grp_vec   <- sub("^[A-D]\\d+_([^_]+)_([MC]).*$", "\\1_\\2", colnames(log_mat))
min_valid <- 2   # of 5 replicates

ok_per_grp <- sapply(unique(grp_vec), function(g)
  rowSums(!is.na(log_mat[, grp_vec == g, drop = FALSE])) >= min_valid)
keep_mv <- rowSums(ok_per_grp) >= 1


mat     <- mat[keep_mv, , drop = FALSE]
log_mat <- log_mat[keep_mv, , drop = FALSE]
df      <- df[df$Protein.Group %in% rownames(mat), ]   # keep df aligned


# design table: one row per sample, 'sample' must match the matrix column names.
# 'group' = line x fraction (8 groups); reused later to build pca_meta.

norm_design <- data.frame(sample = qcols, stringsAsFactors = FALSE)

norm_design$line     <- sub("^[A-D]\\d+_([^_]+)_.*$", "\\1", norm_design$sample)
norm_design$fraction <- sub("^.*_([MC])\\d+.*$",     "\\1", norm_design$sample)
norm_design$rep      <- as.integer(sub("^.*_[MC](\\d+).*$", "\\1", norm_design$sample))
norm_design$group    <- paste(norm_design$line, norm_design$fraction, sep = "_")


stopifnot(nrow(norm_design) == 40,
          all(table(norm_design$group) == 5),
          !any(is.na(norm_design$rep)))

# 4. Annotation sets ----------------------------------------------------------

## MitoCarta 3.0
mitocarta <- read_excel("DIA data/Human.MitoCarta3.0.xls",
                        sheet = "A Human MitoCarta3.0") %>% as.data.frame()

mito_accessions <- mitocarta$UniProt %>%
  as.character() %>%
  strsplit("[|;,]") %>% unlist() %>% trimws() %>%
  sub("-\\d+$", "", .) %>%
  setdiff(c("0", "0.0", "", NA)) %>%
  unique()






## Split a protein group into member accessions; strip isoform suffixes and
## drop contaminant entries, which would otherwise pad the ORA universe.
split_acc <- function(x) {
  strsplit(x, ";") |>
    lapply(function(a) {
      a <- sub("-\\d+$", "", trimws(a))
      a[!grepl("^(Cont_|cRAP-)", a)]
    })
}
# group-aware matching (test every accession in a protein group)

acc_list <- split_acc(rownames(mat))

mito_symbols <- unique(trimws(as.character(mitocarta$Symbol)))
gene_list <- lapply(
  strsplit(ifelse(is.na(df$Genes[match(rownames(mat), df$Protein.Group)]), "",
                  df$Genes[match(rownames(mat), df$Protein.Group)]), ";"),
  trimws)

is_mito <- vapply(acc_list,  function(a) any(a %in% mito_accessions), logical(1)) |
  vapply(gene_list, function(g) any(g %in% mito_symbols),    logical(1))




# 5. PCA plots ----------------------------------------------------------------

run_pca <- function(m, meta) {
  cc  <- m[complete.cases(m), , drop = FALSE]
  pc  <- prcomp(t(cc), scale. = FALSE)          # per-protein centering (default)
  vp  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  list(df = cbind(meta, PC1 = pc$x[,1], PC2 = pc$x[,2], PC3 = pc$x[,3]),
       var = vp, n = nrow(cc))
}

# sample metadata, in the exact order of qcols / matrix columns
pca_meta <- norm_design[, c("sample","line","fraction")]
pca_meta$line <- factor(pca_meta$line, levels = c("I","29","K","30"),
                        labels = c("ISO","CS29","KOLF","CS30"))
pca_meta$fraction <- factor(pca_meta$fraction, levels = c("M","C"),
                            labels = c("Mito","Cyto"))

# Acquisition batch: reps 1-3 and 4-5 were processed in separate rounds.
pca_meta$rep   <- as.integer(sub("^.*_[MC](\\d+)$", "\\1", pca_meta$sample))
pca_meta$batch <- factor(ifelse(pca_meta$rep >= 4, "B2", "B1"))

# Batch must appear in every line, or batch and genotype are confounded.
stopifnot(all(table(pca_meta$line, pca_meta$batch) > 0))



# Enrichment fractionation profile -----------------------------------------



MIN_VALID_ENR <- 2      # observations required in EACH fraction for a ratio
DEPTH_MIN     <- 3000   # proteins per sample below which a fraction is flagged

lines_all <- c("ISO", "CS29", "KOLF", "CS30")

ppm     <- sweep(mat, 2, colSums(mat, na.rm = TRUE), "/") * 1e6
log_ppm <- log2(ppm)

stopifnot(all(abs(colSums(ppm, na.rm = TRUE) - 1e6) < 1))

# -----------------------------------------------------------------------------
# Per-line enrichment table
# -----------------------------------------------------------------------------
enr_one_line <- function(line_name) {
  
  m_cols <- pca_meta$sample[pca_meta$line == line_name & pca_meta$fraction == "Mito"]
  c_cols <- pca_meta$sample[pca_meta$line == line_name & pca_meta$fraction == "Cyto"]
  
  nM <- rowSums(!is.na(log_ppm[, m_cols, drop = FALSE]))
  nC <- rowSums(!is.na(log_ppm[, c_cols, drop = FALSE]))
  
  # Require >= MIN_VALID_ENR in EACH fraction: sparse ratios are MNAR-biased
  # toward apparent enrichment.
  ok <- nM >= MIN_VALID_ENR & nC >= MIN_VALID_ENR
  
  out <- data.frame(
    Protein.Group = rownames(log_ppm)[ok],
    gene    = sub(";.*$", "", df$Genes[match(rownames(log_ppm)[ok], df$Protein.Group)]),
    line    = line_name,
    log2_MC = rowMeans(log_ppm[ok, m_cols, drop = FALSE], na.rm = TRUE) -
      rowMeans(log_ppm[ok, c_cols, drop = FALSE], na.rm = TRUE),
    is_mito = is_mito[ok],
    n_mito_obs = nM[ok], n_cyto_obs = nC[ok],
    stringsAsFactors = FALSE
  )
  
  
  out[is.finite(out$log2_MC), ]
}

enr <- do.call(rbind, lapply(lines_all, enr_one_line))
enr$line  <- factor(enr$line, levels = lines_all)
enr$class <- factor(ifelse(enr$is_mito, "Mitochondrial (MitoCarta3.0)",
                           "Non-mitochondrial"),
                    levels = c("Mitochondrial (MitoCarta3.0)", "Non-mitochondrial"))

# -----------------------------------------------------------------------------
#Statistics
# -----------------------------------------------------------------------------

enr_stats <- do.call(rbind, lapply(lines_all, function(l) {
  d <- enr[enr$line == l, ]
  a <- d$log2_MC[d$is_mito]; b <- d$log2_MC[!d$is_mito]
  w   <- suppressWarnings(wilcox.test(a, b, conf.int = FALSE))
  auc <- unname(w$statistic) / (length(a) * length(b))
  cyto_depth <- max(colSums(!is.na(
    mat[, pca_meta$sample[pca_meta$line == l & pca_meta$fraction == "Cyto"],
        drop = FALSE])))
  data.frame(line = l, n_mito = length(a), n_other = length(b),
             median_mito = median(a), median_other = median(b),
             diff_medians = median(a) - median(b),
             fold = 2^abs(median(a) - median(b)),
             cliffs_delta = 2 * auc - 1, p = w$p.value,
             cyto_depth = cyto_depth, depth_ok = cyto_depth >= DEPTH_MIN)
}))
enr_stats$line <- factor(enr_stats$line, levels = lines_all)
print(enr_stats, digits = 3)

# -----------------------------------------------------------------------------
#Enrichment profile histogram
# -----------------------------------------------------------------------------
pal2     <- c("Mitochondrial (MitoCarta3.0)" = "#1B7A6B",
              "Non-mitochondrial"            = "grey70")
pal_line <- c("Mitochondrial (MitoCarta3.0)" = "#0F4F45",   # KDE curves, darker
              "Non-mitochondrial"            = "#5E6368")

enr_stats$lab <- sprintf("median shift %+.2f (%.1f-fold)\nCliff's delta = %+.2f\np = %s",
                         enr_stats$diff_medians, enr_stats$fold,
                         enr_stats$cliffs_delta, format.pval(enr_stats$p, digits = 2))
enr_stats$flag <- ifelse(enr_stats$depth_ok, "",
                         sprintf("cytosolic fraction under-sampled (%d proteins)",
                                 enr_stats$cyto_depth))

# facet strips carry n, so a sparse panel cannot be mistaken for a solid one
fac_lab <- setNames(sprintf("%s  (n = %s, %s mito)", enr_stats$line,
                            format(enr_stats$n_mito + enr_stats$n_other, big.mark = ","),
                            enr_stats$n_mito), enr_stats$line)

p_enr <- ggplot(enr, aes(log2_MC, fill = class, colour = class)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 position = "identity", alpha = 0.55, colour = NA) +
  geom_density(fill = NA, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  geom_text(data = enr_stats, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.06, vjust = 1.25, size = 2.7, lineheight = 1.05,
            inherit.aes = FALSE) +
  geom_text(data = enr_stats, aes(x = Inf, y = -Inf, label = flag),
            hjust = 1.03, vjust = -0.8, size = 2.4, colour = "firebrick",
            fontface = "italic", inherit.aes = FALSE) +
  scale_fill_manual(values = pal2) +
  scale_colour_manual(values = pal_line, guide = "none") +
  coord_cartesian(xlim = c(-8, 8)) +
  facet_wrap(~ line, nrow = 1, labeller = labeller(line = fac_lab)) +
  labs(x = expression(log[2]*"(Mito / Cyto)"), y = "Density", fill = NULL,
       title = "Mitochondrial enrichment profile",
       subtitle = paste0("ppm-normalised: each sample scaled to 1e6 total signal\n",
                         "proteins require at least ", MIN_VALID_ENR,
                         " observations in each fraction")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 9))

save_fig(p_enr, "enrichment", "enrichment_profile_histogram", w = 13, h = 4)




MIN_OBS_FRAC <- 0.5    # protein must be seen in >= this fraction of samples
# to help define the reference profile

obs      <- rowSums(!is.na(log_mat))
ref_rows <- obs >= MIN_OBS_FRAC * ncol(log_mat)

if (sum(ref_rows) < 200)
  stop("Only ", sum(ref_rows), " proteins observed in >= ",
       100 * MIN_OBS_FRAC, "% of samples - per-sample offsets unreliable. ",
       "Check missingness before normalising.")

ref         <- rowMeans(log_mat[ref_rows, , drop = FALSE], na.rm = TRUE)
size_factor <- apply(log_mat[ref_rows, , drop = FALSE] - ref, 2,
                     median, na.rm = TRUE)
norm_mat    <- sweep(log_mat, 2, size_factor, "-")


## (a) ALL-40 PCA: per sample normalisation -----------------------------------


pa <- run_pca(norm_mat, pca_meta)

p_pca_all <- ggplot(pa$df, aes(PC1, PC2, colour = line, shape = fraction)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 4) +
  scale_colour_npg() +
  labs(title = "PCA - all 40 samples",
       subtitle = paste0("median-of-ratios offsets from ", sum(ref_rows),
                         " proteins; PCA on ", pa$n, " complete cases"),
       x = paste0("PC1 (", pa$var[1], "%)"),
       y = paste0("PC2 (", pa$var[2], "%)"),
       colour = "Line", shape = "Fraction") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))



save_fig(p_pca_all, "pca", "PCA_all40_persample_PC1_PC2")


## (a2) PC2 vs PC3 ------------------------------------------------------------

p_pca_pc23 <- ggplot(pa$df, aes(PC2, PC3, colour = line, shape = fraction)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 4) +
  scale_colour_npg() +
  labs(title = "PCA - all 40 samples, PC2 vs PC3",
       subtitle = "PC1 (fraction) omitted to expose line/disease structure",
       x = paste0("PC2 (", pa$var[2], "%)"),
       y = paste0("PC3 (", pa$var[3], "%)"),
       colour = "Line", shape = "Fraction") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p_pca_pc23, "pca", "PCA_all40_PC2_PC3")


## (b) MITO-ONLY PCA: global median normalisation ----------------------------


mito_cols  <- pca_meta$sample[pca_meta$fraction == "Mito"]
log_mito    <- log_mat[, mito_cols, drop = FALSE]

# ---- ACTIVE: global median normalisation ------------------------------------
# Centre every mito sample on the same median. Assumes most proteins are stable.
mito_med  <- apply(log_mito, 2, median, na.rm = TRUE)
grand_med <- median(mito_med)
mito_glob <- sweep(log_mito, 2, mito_med - grand_med, "-")   # global median norm

# ---- ALTERNATIVE: invariant-subset normalisation ----------------------------
# Anchors on the top 15% most rank-stable proteins. To use, comment out the
# global-median block above and uncomment this.
#
# cc_mito <- complete.cases(log_mito)
# cc_dat  <- log_mito[cc_mito, , drop = FALSE]
#
# rank_mat <- apply(cc_dat, 2, rank)                 # rank within each sample
# rank_sd  <- apply(rank_mat, 1, sd)                 # stability score per protein
#
# frac_keep <- 0.15                                  # top 15% most stable
# n_keep    <- max(200, floor(frac_keep * nrow(cc_dat)))
# invariant <- names(sort(rank_sd))[seq_len(n_keep)]
#
# inv_dat   <- cc_dat[invariant, , drop = FALSE]
# grand_med <- median(inv_dat)
# mito_med  <- apply(inv_dat, 2, median)
# mito_glob <- sweep(log_mito, 2, mito_med - grand_med, "-")   # invariant-subset norm


pm <- run_pca(mito_glob, pca_meta[pca_meta$fraction == "Mito", ])

p_pca_mito <- ggplot(pm$df, aes(PC1, PC2, colour = line)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_colour_npg() +
  labs(title = "PCA - mito fraction only ",
       subtitle = "ISO/CS29 vs KOLF/CS30",
       x = paste0("PC1 (", pm$var[1], "%)"),
       y = paste0("PC2 (", pm$var[2], "%)"),
       colour = "Line") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p_pca_mito, "pca", "PCA_mito_only_PC1_PC2")



# 6.Volcano - two independent line-vs-line models, WITHIN the mito fraction ----

#ISO vs CS29
#KOLF vs CS30 
#Models are run SEPARATELY so the tight isogenic pair and the noisy
#cross-background pair do not share a pooled variance estimate.


mito_norm <- mito_glob                       

## helper: fit one batch-adjusted contrast on a chosen pair of lines,
## return the topTable results plus the fit object.
run_pair <- function(norm_m, meta, line_ref, line_alt, min_valid = 2) {
  cols <- meta$sample[meta$line %in% c(line_ref, line_alt)]
  m    <- norm_m[, cols, drop = FALSE]
  
  grp  <- factor(meta$line[match(cols, meta$sample)],
                 levels = c(line_ref, line_alt))    # ref first
  batch <- droplevels(factor(meta$batch[match(cols, meta$sample)]))
  
  
  
  
  # Keep proteins with >= min_valid observations in EACH line. limma fits
  # row-wise, so partially-observed proteins are tested with less power.
  n_ok <- cbind(ref = rowSums(!is.na(m[, grp == line_ref, drop = FALSE])),
                alt = rowSums(!is.na(m[, grp == line_alt, drop = FALSE])))
  keep <- n_ok[, "ref"] >= min_valid & n_ok[, "alt"] >= min_valid
  m    <- m[keep, , drop = FALSE]
  
  
  # Additive batch term. Balanced across lines, so it is orthogonal to the
  # contrast: removes batch variance without biasing logFC. Costs 1 df.
  if (nlevels(batch) > 1) {
    design <- model.matrix(~ 0 + grp + batch)
  } else {
    warning("Only one batch level for ", line_ref, " vs ", line_alt,
            "; fitting without a batch term.")
    design <- model.matrix(~ 0 + grp)
  }
  colnames(design)[seq_along(levels(grp))] <- levels(grp)
  
  ctr <- makeContrasts(contrasts = paste0(line_alt, "-", line_ref),
                       levels = design)
  
  fit <- eBayes(contrasts.fit(lmFit(m, design), ctr),
                trend = TRUE, robust = TRUE)
  
  res <- topTable(fit, number = Inf, sort.by = "p")
  res$Protein.Group <- rownames(res)
  res$Genes <- df$Genes[match(res$Protein.Group, df$Protein.Group)]
  list(res = res[order(res$P.Value), ], fit = fit)
}

# pca_meta has the labelled factor levels ISO/CS29/KOLF/CS30
mito_meta <- pca_meta[pca_meta$fraction == "Mito", ]
mito_meta$line <- as.character(mito_meta$line)     # drop factor for matching


## Proteins seen in every replicate of one line and none of the other.
## No ratio exists, so limma cannot test them - report them qualitatively.
onoff_hits <- function(norm_m, meta, line_ref, line_alt, ann = df) {
  c_ref <- meta$sample[meta$line == line_ref]
  c_alt <- meta$sample[meta$line == line_alt]
  n_ref <- rowSums(!is.na(norm_m[, c_ref, drop = FALSE]))
  n_alt <- rowSums(!is.na(norm_m[, c_alt, drop = FALSE]))
  
  call <- ifelse(n_ref == length(c_ref) & n_alt == 0, paste("only in", line_ref),
                 ifelse(n_alt == length(c_alt) & n_ref == 0, paste("only in", line_alt), NA))
  
  out <- data.frame(Protein.Group = rownames(norm_m),
                    n_ref = n_ref, n_alt = n_alt, call = call,
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$call), ]
  
  # presence in the OTHER lines too, so "only in X" can be read correctly.
  others <- setdiff(unique(meta$line), c(line_ref, line_alt))
  for (o in others)
    out[[paste0("n_", o)]] <-
    rowSums(!is.na(norm_m[out$Protein.Group,
                          meta$sample[meta$line == o], drop = FALSE]))
  
  out$Genes    <- ann$Genes[match(out$Protein.Group, ann$Protein.Group)]
  out$mean_obs <- rowMeans(norm_m[out$Protein.Group, , drop = FALSE], na.rm = TRUE)
  out[order(-out$mean_obs),
      c("Protein.Group", "Genes", "call", "n_ref", "n_alt",
        paste0("n_", others), "mean_obs")]
}



iso_onoff  <- onoff_hits(mito_norm, mito_meta, "ISO",  "CS29")
kolf_onoff <- onoff_hits(mito_norm, mito_meta, "KOLF", "CS30")


iso_cs29  <- run_pair(mito_norm, mito_meta, "ISO",  "CS29")
kolf_cs30 <- run_pair(mito_norm, mito_meta, "KOLF", "CS30")


## Global multiple-testing correction ------------------------------------------
## BH applied once across both contrasts rather than per contrast: stricter, and
## valid under the positive dependence between contrasts sharing data.
## Per-contrast values kept as adj.P.Val.local.

iso_cs29$res$adj.P.Val.local  <- iso_cs29$res$adj.P.Val
kolf_cs30$res$adj.P.Val.local <- kolf_cs30$res$adj.P.Val

.pooled_q <- p.adjust(c(iso_cs29$res$P.Value, kolf_cs30$res$P.Value),
                      method = "BH")
.n_iso    <- nrow(iso_cs29$res)

iso_cs29$res$adj.P.Val  <- .pooled_q[seq_len(.n_iso)]
kolf_cs30$res$adj.P.Val <- .pooled_q[.n_iso + seq_len(nrow(kolf_cs30$res))]

rm(.pooled_q, .n_iso)








cat("\nISO vs CS29 significant (adj.P<0.05):",
    sum(iso_cs29$res$adj.P.Val < 0.05), "\n")
cat("KOLF vs CS30 significant (adj.P<0.05):",
    sum(kolf_cs30$res$adj.P.Val < 0.05), "\n")


pdf(file.path(fig_dirs[["qc"]], "SA_plots.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plotSA(iso_cs29$fit,  cex = 0.5, main = "ISO vs CS29")
plotSA(kolf_cs30$fit, cex = 0.5, main = "KOLF vs CS30")
dev.off()


# Raw p-value histograms - one per model --------------------------------------

# Expect: a spike near 0 (real hits) on a flat uniform floor.
# a peak near 1, or a U-shape, flags a bad model / wrong normalisation
# a huge spike swamping everything can indicate confounding (KOLF/CS30)



pval_hist <- function(res, ttl) {
  ggplot(res, aes(P.Value)) +
    geom_histogram(bins = 40, fill = "grey40", boundary = 0) +
    theme_classic() +
    labs(title = ttl, x = "raw p-value", y = "proteins")
}



save_fig(pval_hist(iso_cs29$res, "p-value distribution - ISO vs CS29"),
         "qc", "pval_hist_ISO_CS29")
save_fig(pval_hist(kolf_cs30$res, "p-value distribution - KOLF vs CS30"),
         "qc", "pval_hist_KOLF_CS30")


# 7.Volcano plots (one per comparison) ----------------------------------------


## Volcano. y-axis is the BH-adjusted p-value, so height and significance come
## from the same quantity and the dashed line sits at the FDR threshold.
volcano <- function(res, ttl, q = 0.05, n_lab = 30) {
  sig  <- !is.na(res$adj.P.Val) & res$adj.P.Val < q
  top  <- head(res[sig, ], n_lab)
  
  # guard against adj.P.Val == 0 producing an infinite -log10
  res$negLogQ <- -log10(pmax(res$adj.P.Val, .Machine$double.xmin))
  top$negLogQ <- -log10(pmax(top$adj.P.Val, .Machine$double.xmin))
  
  ggplot(res, aes(logFC, negLogQ)) +
    geom_point(aes(colour = sig, alpha = sig), pch = 20, size = 2.5) +
    geom_hline(yintercept = -log10(q), linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = top, aes(label = Genes), size = 2,
                    colour = "grey20", max.overlaps = 20) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "orangered3"),
                        labels = c("ns", sprintf("BH FDR < %.2f", q)), name = NULL) +
    scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 1), guide = "none") +
    labs(x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"(BH-adjusted p)"),
         title = ttl,
         subtitle = sprintf("limma (batch-adjusted), %d proteins tested; %d at global BH FDR < %.2f (%d up, %d down); dashed line = FDR threshold",
                            nrow(res), sum(sig), q,
                            sum(sig & res$logFC > 0), sum(sig & res$logFC < 0))) +
    theme_classic() + theme(legend.position = "top")
}


save_fig(volcano(iso_cs29$res,  "ISO vs CS29"),
         "volcano", "volcano_ISO_vs_CS29", w = 7, h = 6)

save_fig(volcano(kolf_cs30$res, "KOLF vs CS30"),
         "volcano", "volcano_KOLF_vs_CS30", w = 7, h = 6)

# -----------------------------------------------------------------------------
# 7b. Concordance between the two contrasts
# -----------------------------------------------------------------------------

concordance <- function(res_a, res_b, lab_a, lab_b, p_cut = 0.05) {
  a <- res_a[, c("Protein.Group", "Genes", "logFC", "adj.P.Val")]
  b <- res_b[, c("Protein.Group", "logFC", "adj.P.Val")]
  names(a)[3:4] <- paste0(c("logFC_", "adjP_"), lab_a)
  names(b)[2:3] <- paste0(c("logFC_", "adjP_"), lab_b)
  m <- merge(a, b, by = "Protein.Group")
  
  fa <- m[[paste0("logFC_", lab_a)]]; qa <- m[[paste0("adjP_", lab_a)]]
  fb <- m[[paste0("logFC_", lab_b)]]; qb <- m[[paste0("adjP_", lab_b)]]
  
  m$shared     <- qa < p_cut & qb < p_cut & sign(fa) == sign(fb)
  m$rank_score <- -log10(qa) + -log10(qb)
  
  # --- all reporting happens BEFORE m is reordered, so fa/fb/qa/qb still align
  rho <- cor(fa, fb, method = "spearman", use = "complete.obs")
  cat(sprintf("\nConcordance %s vs %s: %d proteins tested in both, Spearman rho = %.3f\n",
              lab_a, lab_b, nrow(m), rho))
  cat(sprintf("  significant in both, same direction: %d (up %d / down %d)\n",
              sum(m$shared), sum(m$shared & fa > 0), sum(m$shared & fa < 0)))
  cat(sprintf("  significant in %s only: %d | %s only: %d | discordant: %d\n",
              lab_a, sum(qa < p_cut & qb >= p_cut),
              lab_b, sum(qb < p_cut & qa >= p_cut),
              sum(qa < p_cut & qb < p_cut & sign(fa) != sign(fb))))
  
  # --- reorder only now; everything below uses m alone
  m <- m[order(-m$rank_score), ]
  
  lab <- head(m[m$shared, ], 20)
  p <- ggplot(m, aes(.data[[paste0("logFC_", lab_a)]],
                     .data[[paste0("logFC_", lab_b)]])) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(aes(colour = shared, alpha = shared), size = 1.8) +
    geom_point(data = m[m$shared, ], colour = "orangered3", size = 2.2) +
    geom_text_repel(data = lab, aes(label = Genes), size = 2.4,
                    colour = "grey15", max.overlaps = 25, min.segment.length = 0) +
    scale_colour_manual(values = c("FALSE" = "grey75", "TRUE" = "orangered3")) +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
    labs(x = paste0("log2 FC  ", lab_a), y = paste0("log2 FC  ", lab_b),
         title = "Concordance between the two contrasts",
         subtitle = sprintf("Spearman rho = %.2f; %d significant in both, same direction",
                            rho, sum(m$shared))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          plot.title = element_text(face = "bold"))
  
  save_fig(p, "concordance", "logFC_concordance", w = 6.5, h = 6)
  invisible(list(table = m, rho = rho))
}

con <- concordance(iso_cs29$res, kolf_cs30$res, "CS29_vs_ISO", "CS30_vs_KOLF")


# QC for each model -----------------------------------------------------------




sig_table <- function(res, als_line, ctrl_line, p_cut = 0.05) {
  keep_cols <- c("comparison","Protein.Group","Genes","logFC","fold_change",
                 "direction","AveExpr","P.Value","adj.P.Val")
  
  out <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < p_cut, , drop = FALSE]
  
  # Nothing passes: return an empty table with the right columns rather than
  # erroring, since scalar assignment to a 0-row data frame fails.
  if (nrow(out) == 0) {
    empty <- data.frame(comparison   = character(0),
                        Protein.Group = character(0),
                        Genes         = character(0),
                        logFC         = numeric(0),
                        fold_change   = numeric(0),
                        direction     = character(0),
                        AveExpr       = numeric(0),
                        P.Value       = numeric(0),
                        adj.P.Val     = numeric(0),
                        stringsAsFactors = FALSE)
    return(empty)
  }
  
  out <- out[order(out$adj.P.Val), , drop = FALSE]
  out$comparison  <- paste(als_line, "vs", ctrl_line)
  out$direction   <- ifelse(out$logFC > 0,
                            paste("Up in", als_line),
                            paste("Down in", als_line))
  out$fold_change <- round(2^abs(out$logFC), 2)   # linear fold change
  out[, keep_cols, drop = FALSE]
}

iso_sig  <- sig_table(iso_cs29$res,  als_line = "CS29", ctrl_line = "ISO")
kolf_sig <- sig_table(kolf_cs30$res, als_line = "CS30", ctrl_line = "KOLF")

cat("ISO vs CS29  significant:", nrow(iso_sig),  "\n")
cat("KOLF vs CS30 significant:", nrow(kolf_sig), "\n")

xlsx_path <- file.path(fig_dirs[["tables"]], "DE_significant_mito.xlsx")


writexl::write_xlsx(list("ISO_vs_CS29"        = iso_sig,
                         "KOLF_vs_CS30"       = kolf_sig,
                         "ISO_vs_CS29_onoff"  = iso_onoff,
                         "KOLF_vs_CS30_onoff" = kolf_onoff,
                         "concordance_shared" = con$table[con$table$shared, ],
                         "concordance_all"    = con$table),
                    xlsx_path)


# 8. GO over-representation analysis (ORA) on the DE results ------------------


# take a limma result table, run enrichGO on one direction and one ontology
run_go <- function(res, direction = c("up","down"), ont = "BP", p_cut = 0.05) {
  direction <- match.arg(direction)
  
  # ORA unit = protein group, matching limma. One accession per group, or a
  # single observation would contribute multiple counts to the test.
  rep_acc <- vapply(split_acc(res$Protein.Group),
                    function(a) if (length(a)) a[1] else NA_character_,
                    character(1))
  
  sig <- if (direction == "up") res$adj.P.Val < p_cut & res$logFC > 0
  else                   res$adj.P.Val < p_cut & res$logFC < 0
  
  genes    <- unique(na.omit(rep_acc[sig]))
  universe <- unique(na.omit(rep_acc))
  
  if (length(genes) < 10) {
    warning("Fewer than 10 significant proteins - ORA will be underpowered.")
    return(NULL)
  }
  
  enrichGO(gene         = genes,
           universe     = universe,
           OrgDb        = org.Hs.eg.db,
           keyType      = "UNIPROT",
           ont          = ont,
           pvalueCutoff = 0.05,
           readable     = TRUE)
}

## run all three ontologies for one comparison/direction
run_go_all <- function(res, direction, onts = c("BP","MF","CC")) {
  out <- lapply(onts, function(o) {
    e <- run_go(res, direction, o)
    e
  })
  setNames(out, onts)
}

go_iso_up   <- run_go_all(iso_cs29$res,  "up")
go_iso_down <- run_go_all(iso_cs29$res,  "down")

go_kolf_up   <- run_go_all(kolf_cs30$res, "up")
go_kolf_down <- run_go_all(kolf_cs30$res, "down")

# 8.1 Reduce GO redundancy with semantic similarity ---------------------------


gd_list <- list(BP = godata('org.Hs.eg.db', ont = "BP"),
                MF = godata('org.Hs.eg.db', ont = "MF"),
                CC = godata('org.Hs.eg.db', ont = "CC"))

simplify_go <- function(ego, ont, cutoff = 0.7) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  sim <- pairwise_termsim(ego, method = "Wang", semData = gd_list[[ont]])
  list(sim = sim,
       simplified = simplify(sim, cutoff = cutoff, by = "p.adjust", select_fun = min))
}

## simplify every ontology in a run_go_all() result
simplify_all <- function(go_list, cutoff = 0.7)
  setNames(lapply(names(go_list),
                  function(o) simplify_go(go_list[[o]], o, cutoff)),
           names(go_list))

go_iso_up_s    <- simplify_all(go_iso_up)
go_iso_down_s  <- simplify_all(go_iso_down)
go_kolf_up_s   <- simplify_all(go_kolf_up)
go_kolf_down_s <- simplify_all(go_kolf_down)


# 8.2 Faceted dotplot - top N terms PER ONTOLOGY -----------

GO_TOP_N <- 12   # terms shown per ontology; display only. Change if you wanna see top 15 etc proteins

## collect the simplified results into one tidy table
go_table <- function(s_list, top_n = GO_TOP_N) {
  d <- do.call(rbind, lapply(names(s_list), function(o) {
    s <- s_list[[o]]
    if (is.null(s)) return(NULL)
    x <- as.data.frame(s$simplified)
    if (!nrow(x)) return(NULL)
    x$ONTOLOGY <- o
    head(x[order(x$p.adjust), ], top_n)
  }))
  if (is.null(d) || !nrow(d)) return(NULL)
  d$ONTOLOGY <- factor(d$ONTOLOGY, levels = c("BP","MF","CC"))
  d[order(d$ONTOLOGY, d$p.adjust), ]
}

go_dot <- function(d, ttl) {
  if (is.null(d) || !nrow(d)) {
    message("No enriched terms for: ", ttl); return(invisible(NULL))
  }
  d$label <- make.unique(d$Description)          
  d$label <- factor(d$label, levels = rev(d$label))
  ggplot(d, aes(Count, label, colour = p.adjust, size = Count)) +
    geom_point() +
    scale_colour_gradient(low = "orangered3", high = "grey70",
                          name = "p.adjust") +
    scale_size_continuous(guide = "none") +
    scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Proteins in term", y = NULL, title = ttl,
         subtitle = sprintf("top %d simplified terms per ontology", GO_TOP_N)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey95", colour = NA))
}

d_iso_up    <- go_table(go_iso_up_s)
d_iso_down  <- go_table(go_iso_down_s)
d_kolf_up   <- go_table(go_kolf_up_s)
d_kolf_down <- go_table(go_kolf_down_s)



save_fig(go_dot(d_iso_up,    "GO ORA - up in CS29 vs ISO"),
         "go", "GO_ISO_CS29_up_simplified",    w = 8, h = 10)
save_fig(go_dot(d_iso_down,  "GO ORA - down in CS29 vs ISO"),
         "go", "GO_ISO_CS29_down_simplified",  w = 8, h = 10)
save_fig(go_dot(d_kolf_up,   "GO ORA - up in CS30 vs KOLF"),
         "go", "GO_KOLF_CS30_up_simplified",   w = 8, h = 10)
save_fig(go_dot(d_kolf_down, "GO ORA - down in CS30 vs KOLF"),
         "go", "GO_KOLF_CS30_down_simplified", w = 8, h = 10)






writeLines(capture.output(sessionInfo()),
           file.path(fig_dirs[["tables"]], "sessionInfo.txt"))