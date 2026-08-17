#!/usr/bin/env Rscript
# =============================================================================
# enrichment_distribution.R
# -----------------------------------------------------------------------------
# Answers three questions the deck currently asserts without showing:
#   (1) How strong are the binding events, per strain?  mean / median / range
#   (2) Does enrichment differ by ChExMix SUBTYPE?
#   (3) Does enrichment track REPRODUCIBILITY (replication code)?
#
# Run from the repo root:
#   conda activate multiomics-r      # or the figures env
#   Rscript workflow/scripts/enrichment_distribution.R
#
# Outputs
#   results/07_figures/enrichment_summary.tsv        <- the tables
#   results/07_figures/enrichment_distribution.pdf   <- 3-panel figure
#   results/07_figures/enrichment_distribution.png
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(forcats)
})

## Paths come from the rule, not from hard-coded strings: PEAKS moved under a
## per-caller subdirectory (03_peaks/chexmix/) when a second caller was added,
## and a hard-coded path would have silently read a stale tree.
events_files <- unlist(snakemake@input[["events"]])
repcode_files <- unlist(snakemake@input[["repcodes"]])
STRAINS <- unlist(snakemake@params[["strains"]])
OUTDIR <- snakemake@params[["outdir"]]
logp <- snakemake@log[[1]]
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)

lg <- file(logp, open = "wt")
sink(lg, type = "message")
on.exit({sink(type = "message"); close(lg)}, add = TRUE)

names(events_files) <- STRAINS
names(repcode_files) <- STRAINS

# theme matching the project's other figures ---------------------------------
INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"
MUTE <- "#6B7A8A"; SAND <- "#C9743F"; GREY <- "#B4BEC7"
theme_perr <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#EDF1F3", linewidth = .3),
          axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", hjust = 0),
          plot.subtitle = element_text(colour = MUTE, size = 8.5, hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(8, 12, 8, 8))
}

# ---------------------------------------------------------------------------
# ChExMix .events: header lines start with '#'. The column names live on the
# line beginning '#Point'. The fold-change column is named <strain>_log2Fold,
# so it is located BY NAME rather than by position.
# ---------------------------------------------------------------------------
read_events <- function(strain) {
  path <- events_files[[strain]]
  ln   <- readLines(path, warn = FALSE)

  hdr_line <- grep("^#Point", ln, value = TRUE)
  if (!length(hdr_line)) stop("no '#Point' header line in ", path)
  cols <- strsplit(sub("^#", "", hdr_line[1]), "\t")[[1]]

  dat <- ln[!grepl("^#", ln) & nzchar(ln)]
  if (!length(dat)) stop("no data rows in ", path)

  m  <- do.call(rbind, strsplit(dat, "\t", fixed = TRUE))
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- cols[seq_len(ncol(df))]

  lfc_col <- grep("_log2Fold$", names(df), value = TRUE)[1]
  if (is.na(lfc_col)) stop("no *_log2Fold column found in ", path)

  tibble::tibble(
    strain  = strain,
    point   = trimws(df[["Point"]]),
    log2FC  = suppressWarnings(as.numeric(df[[lfc_col]])),
    subtype = if ("SubtypeName" %in% names(df)) trimws(df[["SubtypeName"]]) else "NA"
  ) |> dplyr::filter(is.finite(log2FC))
}

# ---------------------------------------------------------------------------
# ChExMix replicationcodes.table: column 1 = event Point, column 2 = a code
# indicating which replicates called the event. The code alphabet is not fixed
# across ChExMix versions, so we do NOT hard-code it: we report whatever codes
# are present, and derive an n_reps count from the number of non-'X' characters
# when the code looks like a per-replicate string.
# ---------------------------------------------------------------------------
read_repcodes <- function(strain) {
  path <- repcode_files[[strain]]
  if (!file.exists(path)) {
    warning("missing ", path); return(NULL)
  }
  ln  <- readLines(path, warn = FALSE)
  dat <- ln[!grepl("^#", ln) & nzchar(ln)]
  m   <- do.call(rbind, strsplit(dat, "\t", fixed = TRUE))
  if (is.null(m) || ncol(m) < 2) return(NULL)

  tibble::tibble(
    point = trimws(m[, 1]),
    code  = trimws(m[, 2])
  )
}

# ---------------------------------------------------------------------------
# assemble
# ---------------------------------------------------------------------------
ev <- dplyr::bind_rows(lapply(STRAINS, read_events))
rc <- dplyr::bind_rows(lapply(STRAINS, function(s) {
  x <- read_repcodes(s); if (is.null(x)) return(NULL); dplyr::mutate(x, strain = s)
}))

if (!is.null(rc) && nrow(rc)) {
  ev <- dplyr::left_join(ev, rc, by = c("strain", "point"))
} else {
  ev$code <- NA_character_
}

# Derive a replicate count from the code where the alphabet allows it.
# Codes are reported verbatim as well, so nothing is hidden by this guess.
ev <- ev |>
  dplyr::mutate(
    n_reps = dplyr::case_when(
      is.na(code)   ~ NA_integer_,
      code == "X"   ~ NA_integer_,
      TRUE          ~ as.integer(nchar(gsub("[^1A-Za-z]", "", code)))
    ),
    rep_label = ifelse(is.na(code), "unassigned", code)
  )

message("distinct replication codes observed: ",
        paste(sort(unique(na.omit(ev$code))), collapse = ", "))

# ---------------------------------------------------------------------------
# summary tables
# ---------------------------------------------------------------------------
describe <- function(d, ...) {
  d |> dplyr::group_by(...) |>
    dplyr::summarise(
      n       = dplyr::n(),
      mean    = round(mean(log2FC), 3),
      median  = round(median(log2FC), 3),
      sd      = round(sd(log2FC), 3),
      min     = round(min(log2FC), 3),
      q25     = round(quantile(log2FC, .25), 3),
      q75     = round(quantile(log2FC, .75), 3),
      max     = round(max(log2FC), 3),
      fold_median = round(2^median(log2FC), 2),
      .groups = "drop"
    )
}

by_strain  <- describe(ev, strain) |> dplyr::mutate(stratum = "strain", level = strain)
by_subtype <- describe(ev, strain, subtype) |> dplyr::mutate(stratum = "subtype", level = subtype)
by_reps    <- describe(dplyr::filter(ev, !is.na(code)), strain, rep_label) |>
              dplyr::mutate(stratum = "replication_code", level = rep_label)

summary_tbl <- dplyr::bind_rows(
  dplyr::select(by_strain,  stratum, strain, level, n, mean, median, sd, min, q25, q75, max, fold_median),
  dplyr::select(by_subtype, stratum, strain, level, n, mean, median, sd, min, q25, q75, max, fold_median),
  dplyr::select(by_reps,    stratum, strain, level, n, mean, median, sd, min, q25, q75, max, fold_median)
)

readr_out <- file.path(OUTDIR, "enrichment_summary.tsv")
write.table(summary_tbl, readr_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", readr_out)
print(as.data.frame(summary_tbl), row.names = FALSE)

# ---------------------------------------------------------------------------
# PANEL 1 — distribution per strain, with mean and median marked
# ---------------------------------------------------------------------------
stats1 <- describe(ev, strain)
minfold_line <- log2(1.5)   # the --minfold 1.5 floor, in log2 units

p1 <- ggplot(ev, aes(log2FC, fill = strain, colour = strain)) +
  geom_density(alpha = .28, linewidth = .6, adjust = .9) +
  geom_vline(xintercept = minfold_line, colour = GREY, linetype = "dotted", linewidth = .5) +
  geom_vline(data = stats1, aes(xintercept = median, colour = strain),
             linetype = "solid", linewidth = .7) +
  geom_vline(data = stats1, aes(xintercept = mean, colour = strain),
             linetype = "dashed", linewidth = .55) +
  annotate("text", x = minfold_line, y = Inf, label = "  --minfold 1.5 floor",
           hjust = 0, vjust = 1.6, size = 2.8, colour = MUTE) +
  scale_fill_manual(values = c(A909 = RUST, CJB111 = TEAL)) +
  scale_colour_manual(values = c(A909 = RUST, CJB111 = TEAL)) +
  labs(title = "Enrichment distribution of called binding events",
       subtitle = "solid line = median   ·   dashed = mean   ·   dotted = the caller's minimum-fold floor",
       x = "ChExMix log2 fold enrichment over scaled control", y = "density") +
  theme_perr()

# ---------------------------------------------------------------------------
# PANEL 2 — by ChExMix subtype (A909 only has >1 subtype)
# ---------------------------------------------------------------------------
p2 <- ggplot(ev, aes(x = fct_reorder(subtype, log2FC, .fun = median),
                     y = log2FC, fill = strain)) +
  geom_boxplot(outlier.size = .8, alpha = .55, width = .6, linewidth = .4) +
  geom_jitter(width = .12, height = 0, size = .9, alpha = .55, colour = INK) +
  geom_hline(yintercept = minfold_line, colour = GREY, linetype = "dotted") +
  facet_wrap(~ strain, scales = "free_x") +
  scale_fill_manual(values = c(A909 = RUST, CJB111 = TEAL), guide = "none") +
  labs(title = "Enrichment by ChExMix binding subtype",
       subtitle = "subtype identity comes from the tag-distribution SHAPE, not from the label",
       x = NULL, y = "log2 fold enrichment") +
  theme_perr() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ---------------------------------------------------------------------------
# PANEL 3 — by replication code (reproducibility)
# ---------------------------------------------------------------------------
ev_rc <- dplyr::filter(ev, !is.na(code))
if (nrow(ev_rc)) {
  p3 <- ggplot(ev_rc, aes(x = fct_reorder(rep_label, log2FC, .fun = median),
                          y = log2FC, fill = strain)) +
    geom_boxplot(outlier.size = .8, alpha = .55, width = .6, linewidth = .4) +
    geom_jitter(width = .12, height = 0, size = .9, alpha = .55, colour = INK) +
    geom_hline(yintercept = minfold_line, colour = GREY, linetype = "dotted") +
    facet_wrap(~ strain, scales = "free_x") +
    scale_fill_manual(values = c(A909 = RUST, CJB111 = TEAL), guide = "none") +
    labs(title = "Enrichment by replication code",
         subtitle = "if enrichment rises with reproducibility, 'called in all replicates' is a principled filter",
         x = "ChExMix replication code", y = "log2 fold enrichment") +
    theme_perr() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
} else {
  p3 <- ggplot() + theme_void() +
    annotate("text", 0, 0, label = "no replication codes joined", colour = MUTE)
}

fig <- (p1 / (p2 | p3)) +
  plot_annotation(
    title = "How strong are the binding events, really?",
    subtitle = "Full distribution, and the two stratifications that decide how much weight each event carries",
    theme = theme(plot.title = element_text(face = "bold", size = 15, colour = INK),
                  plot.subtitle = element_text(colour = MUTE, size = 10)))

ggsave(file.path(OUTDIR, "enrichment_distribution.pdf"), fig,
       width = 12, height = 9.5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "enrichment_distribution.png"), fig,
       width = 12, height = 9.5, dpi = 300)
message("wrote ", file.path(OUTDIR, "enrichment_distribution.pdf"))