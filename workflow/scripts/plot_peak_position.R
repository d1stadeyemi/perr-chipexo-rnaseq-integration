## =============================================================================
## plot_peak_position.R -- where do peaks sit relative to genes?
## =============================================================================
## Replaces a three-category table (inside / upstream / operon) with the
## underlying continuous distribution. The categories are thresholds we chose;
## the distribution is a property of the data, and it answers directly the
## question the table invites: "why are most peaks inside coding sequence?"
##
## Bacterial genomes are ~88% coding, so a peak landing in a CDS is the null
## expectation, not an anomaly. Panel B makes that explicit by comparing the
## observed fraction against the coding fraction of each genome.
## =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(scales)
})

assign_files <- unlist(snakemake@input[["assign"]])
gff_files    <- unlist(snakemake@input[["gff"]])
STRAINS      <- unlist(snakemake@params[["strains"]])
MAX_UP       <- as.integer(snakemake@params[["max_upstream"]])
outdir       <- snakemake@params[["outdir"]]
logp         <- snakemake@log[[1]]

names(assign_files) <- STRAINS; names(gff_files) <- STRAINS
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt"); sink(lg, type = "message")
on.exit({sink(type = "message"); close(lg)}, add = TRUE)

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"
MUTE <- "#6B7A8A"; GREY <- "#B4BEC7"; PALE <- "#EDF1F3"; SAND <- "#C9743F"
SCOL <- c(A909 = RUST, CJB111 = TEAL)

theme_perr <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = PALE, linewidth = .3),
          axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", hjust = 0, size = 12),
          plot.subtitle = element_text(colour = MUTE, size = 8.5, hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(8, 12, 8, 8))
}

## ---------------------------------------------------------------------------
## coding fraction per genome -- the null expectation for "inside"
## ---------------------------------------------------------------------------
coding_fraction <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!startsWith(ln, "#")]
  f  <- strsplit(ln, "\t", fixed = TRUE)
  f  <- f[lengths(f) >= 9]
  ty <- vapply(f, `[`, "", 3)
  st <- as.numeric(vapply(f, `[`, "", 4)); en <- as.numeric(vapply(f, `[`, "", 5))
  glen <- max(en[ty == "region"], na.rm = TRUE)
  keep <- ty == "gene"
  ## merge overlapping genes so the fraction cannot exceed 1
  iv <- data.frame(s = st[keep], e = en[keep])
  iv <- iv[order(iv$s), ]
  tot <- 0; cs <- iv$s[1]; ce <- iv$e[1]
  for (i in seq_len(nrow(iv))[-1]) {
    if (iv$s[i] <= ce) { ce <- max(ce, iv$e[i]) } else { tot <- tot + (ce - cs); cs <- iv$s[i]; ce <- iv$e[i] }
  }
  tot <- tot + (ce - cs)
  tot / glen
}

## ---------------------------------------------------------------------------
## one row per PEAK, not per peak-gene pair: a peak assigned to three genes
## would otherwise be counted three times
## ---------------------------------------------------------------------------
dat <- bind_rows(lapply(STRAINS, function(s) {
  d <- suppressMessages(read_tsv(assign_files[[s]], show_col_types = FALSE))
  d |>
    mutate(dist = suppressWarnings(as.numeric(dist_to_TSS))) |>
    filter(!is.na(dist), relationship %in% c("inside", "upstream")) |>
    group_by(peak_pos) |>
    ## nearest assignment defines the peak's position
    slice_min(abs(dist), n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(strain = s, dist, relationship)
}))

cf <- tibble(strain = STRAINS,
             coding = vapply(STRAINS, function(s) coding_fraction(gff_files[[s]]), 0))
for (k in seq_len(nrow(cf)))
  message(sprintf("%s: coding fraction of genome = %.1f%%", cf$strain[k], 100 * cf$coding[k]))

obs <- dat |> group_by(strain) |>
  summarise(n = n(), inside = sum(relationship == "inside"),
            pct_inside = 100 * inside / n, .groups = "drop") |>
  left_join(cf, by = "strain") |> mutate(expected = 100 * coding)
for (k in seq_len(nrow(obs)))
  message(sprintf("%s: %d peaks, %.0f%% inside a gene (expected ~%.0f%% by chance)",
                  obs$strain[k], obs$n[k], obs$pct_inside[k], obs$expected[k]))

## ---------------------------------------------------------------------------
## PANEL A -- the distribution
## ---------------------------------------------------------------------------
pA <- ggplot(dat, aes(dist, fill = strain)) +
  annotate("rect", xmin = -MAX_UP, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = PALE, alpha = .9) +
  geom_histogram(binwidth = 50, colour = "white", linewidth = .2, alpha = .85) +
  geom_vline(xintercept = 0, colour = INK, linewidth = .5) +
  ## Shared y-axis so the two strains are directly comparable. A909 has ~4x as
  ## many peaks, so CJB111's bars are correspondingly short -- that difference
  ## is real and should not be hidden by rescaling each panel to its own maximum.
  facet_wrap(~ strain, ncol = 1) +
  scale_fill_manual(values = SCOL, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(8)) +
  labs(title = "A   Peak position relative to the nearest transcription start",
       subtitle = paste0("shaded = the ", MAX_UP,
                         " bp promoter window;  0 = TSS;  positive = inside the gene"),
       x = "distance from TSS (bp)", y = "peaks") +
  theme_perr()

## ---------------------------------------------------------------------------
## PANEL B -- observed vs expected, so "inside" is not read as anomalous
## ---------------------------------------------------------------------------
bar <- obs |>
  select(strain, observed = pct_inside, expected) |>
  pivot_longer(-strain, names_to = "what", values_to = "pct") |>
  mutate(what = factor(what, levels = c("expected", "observed"),
                       labels = c("coding fraction of the genome",
                                  "peaks inside a gene")))
pB <- ggplot(bar, aes(strain, pct, fill = what)) +
  geom_col(position = position_dodge(width = .72), width = .58) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = .72), vjust = -0.4,
            size = 3, colour = INK) +
  scale_fill_manual(values = c("coding fraction of the genome" = GREY,
                               "peaks inside a gene" = SAND)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .12)),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "B   Peaks in coding sequence are the expectation",
       subtitle = "bacterial genomes are densely coding, so an intragenic peak is not anomalous",
       x = NULL, y = NULL) +
  theme_perr()

fig <- (pA | pB) + plot_layout(widths = c(1.5, 1)) +
  plot_annotation(
    title = "Where PerR binds relative to genes",
    subtitle = "One point per peak: a peak assigned to several genes is counted once, at its nearest assignment",
    theme = theme(plot.title = element_text(face = "bold", size = 15, colour = INK),
                  plot.subtitle = element_text(colour = MUTE, size = 10)))

ggsave(file.path(outdir, "peak_position.pdf"), fig, width = 12, height = 6.5)
ggsave(file.path(outdir, "peak_position.png"), fig, width = 12, height = 6.5, dpi = 300)
message("wrote peak_position.pdf")