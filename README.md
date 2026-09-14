# Mitochondrial proteomics of isogenic ALS lines

Label-free DIA-MS analysis of mitochondrial and cytosolic fractions from two
independent isogenic pairs of iPSC-derived lines: **CS29 vs ISO** and
**CS30 vs KOLF** (ALS vs gene-corrected control).

## Data

40 samples — 4 lines × 2 fractions (Mito, Cyto) × 5 biological replicates —
quantified with DIA-NN and supplied as a protein-group matrix. Replicates 1–3
and 4–5 were acquired in separate batches.

## Pipeline (`DIA_LFQ_R.R`)

1. **Filtering** — contaminant removal; proteins retained with ≥2 non-missing
   values in at least one group.
2. **Normalisation** — per-sample median-of-ratios against a reference set of
   proteins observed in ≥50% of samples.
3. **QC** — PCA, per-sample intensity distributions, p-value histograms,
   limma mean–variance (SA) plots.
4. **Fraction enrichment** — mitochondrial enrichment assessed against
   MitoCarta 3.0 to confirm fractionation quality.
5. **Differential abundance** — `limma` within the mitochondrial fraction.
   Design `~ 0 + genotype + batch`; batch is balanced across genotypes, so it
   is orthogonal to the contrast and removes technical variance without biasing
   fold changes. Empirical Bayes moderation (`trend`, `robust`).
6. **Multiple testing** — Benjamini–Hochberg applied once across both
   contrasts; volcano plots show BH-adjusted p-values on the y-axis.
7. **Enrichment** — GO over-representation analysis (`clusterProfiler`) on
   significant proteins, split by direction.

## Sensitivity analysis

`quickcheck.R` repeats the pipeline using replicates 1–3 only (batch excluded
by design rather than modelled) with invariant-subset normalisation.

## Requirements

R ≥ 4.2 — `limma`, `tidyverse`, `clusterProfiler`, `org.Hs.eg.db`,
`enrichplot`, `readxl`, `writexl`, `ggrepel`, `ggsci`.
