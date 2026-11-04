## =============================================================================
## plot_motif_evidence.R -- the sequence evidence for the GBS PerR motif
## =============================================================================
## Writes FOUR SEPARATE figures (no combined panel, no titles), sized for
## publication. Panel letters and captions are added in the figure layout.
##
##   motif_logo        what does the motif look like      -> sequence logo
##   motif_position    is it at the binding event          -> offset distribution,
##                                                            GBS motif vs GAS box
##   motif_gbs_vs_gas  like-for-like event counts          -> GBS motif vs GAS box
##   motif_occupancy   does motif presence predict binding -> genome-wide vs bound
##
## Each is written as .pdf (vector, for the paper) and .png (300 dpi, for slides).
## =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(ggseqlogo); library(scales)
})

class_files <- unlist(snakemake@input[["classification"]])
pwm_file    <- snakemake@input[["pwm"]]
fimo_gbs    <- unlist(snakemake@input[["fimo_gbs"]])
fimo_gas    <- unlist(snakemake@input[["fimo_gas"]])
STRAINS     <- unlist(snakemake@params[["strains"]])
HALFWIN     <- as.integer(snakemake@params[["halfwindow"]])
outdir      <- snakemake@params[["outdir"]]
logp        <- snakemake@log[[1]]

names(class_files) <- STRAINS; names(fimo_gbs) <- STRAINS; names(fimo_gas) <- STRAINS
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt"); sink(lg, type = "message")

## ===========================================================================
## FONT AND SIZE CONTROLS -- change these to resize everything
## ===========================================================================
BASE_SIZE  <- 22    # tick labels, legend text
TITLE_SIZE <- 24    # axis titles and strain (facet) labels
BAR_LABEL  <- 7     # numbers printed on bars (mm; 7 is about 20 pt)
FIG_W      <- 8     # default width of each saved figure (inches)
FIG_H      <- 6.5   # default height of each saved figure (inches)
## ===========================================================================

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"
GREY <- "#B4BEC7"; BLUE <- "#606064"; PALE <- "#EDF1F3"
MCOL <- c("GBS motif" = RUST, "GAS Per box" = BLUE)

theme_pub <- function() {
  theme_minimal(base_size = BASE_SIZE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = PALE, linewidth = .4),
          text        = element_text(colour = INK, face = "bold"),
          axis.text   = element_text(colour = INK, face = "bold", size = BASE_SIZE),
          axis.title  = element_text(colour = INK, face = "bold", size = TITLE_SIZE),
          strip.text  = element_text(colour = INK, face = "bold", size = TITLE_SIZE),
          legend.text = element_text(colour = INK, face = "bold", size = BASE_SIZE),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          plot.margin = margin(12, 16, 12, 12))
}

save_fig <- function(p, stem, w = FIG_W, h = FIG_H) {
  ggsave(file.path(outdir, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  message("wrote ", stem)
}

## ---------------------------------------------------------------------------
## 1. sequence logo, read straight from the MEME PWM that FIMO scanned
## ---------------------------------------------------------------------------
read_meme_pwm <- function(path) {
  ln <- readLines(path, warn = FALSE)
  i  <- grep("letter-probability matrix", ln)[1]
  if (is.na(i)) stop("no letter-probability matrix in ", path)
  w  <- as.integer(sub(".*w=\\s*([0-9]+).*", "\\1", ln[i]))
  m  <- do.call(rbind, lapply(ln[(i + 1):(i + w)], function(x)
        as.numeric(strsplit(trimws(x), "\\s+")[[1]])))
  t(structure(m, dimnames = list(NULL, c("A", "C", "G", "T"))))
}
pwm <- read_meme_pwm(pwm_file)
message("PWM width ", ncol(pwm), "  majority consensus ",
        paste(c("A","C","G","T")[apply(pwm, 2, which.max)], collapse = ""))

p_logo <- ggseqlogo(pwm, method = "bits") +
  labs(x = "Position", y = "Bits") +
  theme_pub() + theme(legend.position = "none")
save_fig(p_logo, "motif_logo", w = 8, h = 4.5)

## ---------------------------------------------------------------------------
## 2. position relative to the binding event: GBS motif vs GAS Per box
## For each event, the best-scoring (lowest p) hit of each motif inside the
## event's discovery window is used. The same rule is applied to both motifs,
## so the two distributions are directly comparable.
## ---------------------------------------------------------------------------
read_fimo_hits <- function(path) {
  f <- suppressMessages(read_tsv(path, comment = "#", show_col_types = FALSE))
  f <- f[!is.na(f$start), ]
  tibble(chrom  = sub("^chr", "", f$sequence_name),
         centre = (f$start + f$stop) / 2,
         p      = f$`p-value`)
}

best_offsets <- function(cls, hits, label, strain) {
  cls <- cls |>
    mutate(chrom = sub("^chr", "", chrom),
           point = suppressWarnings(as.numeric(sub(".*_", "", peak))))
  out <- lapply(seq_len(nrow(cls)), function(k) {
    h <- hits[hits$chrom == cls$chrom[k] &
              hits$centre >= cls$window_start[k] &
              hits$centre <= cls$window_end[k], ]
    if (nrow(h) == 0) return(NULL)
    tibble(offset = h$centre[which.min(h$p)] - cls$point[k])
  })
  bind_rows(out) |> mutate(motif = label, strain = strain)
}

off <- bind_rows(lapply(STRAINS, function(s) {
  cls <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  bind_rows(best_offsets(cls, read_fimo_hits(fimo_gbs[[s]]), "GBS motif",   s),
            best_offsets(cls, read_fimo_hits(fimo_gas[[s]]), "GAS Per box", s))
})) |>
  mutate(motif = factor(motif, levels = c("GBS motif", "GAS Per box")))

med <- off |> filter(motif == "GBS motif") |>
  group_by(strain) |> summarise(m = median(offset), n = n(), .groups = "drop")
for (k in seq_len(nrow(med)))
  message(sprintf("%s: GBS hits in window %d, median offset %+d bp",
                  med$strain[k], med$n[k], as.integer(med$m[k])))
gas_n <- off |> filter(motif == "GAS Per box") |> count(strain)
message("GAS hits in window: ",
        paste(sprintf("%s %d", gas_n$strain, gas_n$n), collapse = " | "))

for (s in STRAINS) {
  d  <- off |> filter(strain == s)
  ms <- med |> filter(strain == s)
  p_pos <- ggplot(d, aes(offset, fill = motif)) +
    annotate("rect", xmin = -HALFWIN, xmax = HALFWIN, ymin = -Inf, ymax = Inf,
             fill = PALE, alpha = .9) +
    geom_histogram(binwidth = 10, colour = "white", linewidth = .2,
                   position = "identity", alpha = .9) +
    geom_vline(xintercept = 0, colour = INK, linewidth = .6) +
    geom_vline(data = ms, aes(xintercept = m), colour = RUST,
               linetype = "dashed", linewidth = .9) +
    scale_fill_manual(values = MCOL) +
    scale_x_continuous(breaks = pretty_breaks(5)) +
    scale_y_continuous(breaks = function(l) unique(floor(pretty(c(0, l[2]), 4)))) +
    labs(x = "Distance from binding event centre (bp)",
         y = "Motif occurrences") +
    theme_pub()
  save_fig(p_pos, paste0("motif_position_", s), w = 9, h = 5.5)
}

## ---------------------------------------------------------------------------
## 3. events containing each motif, judged by the identical rule
## ---------------------------------------------------------------------------
cmp <- bind_rows(lapply(STRAINS, function(s) {
  d <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  tibble(strain = s, n_peaks = nrow(d),
         `GBS motif`   = sum(d$gbs_perbox_in_window == TRUE, na.rm = TRUE),
         `GAS Per box` = sum(d$perbox_in_window     == TRUE, na.rm = TRUE))
})) |>
  pivot_longer(c(`GBS motif`, `GAS Per box`), names_to = "motif", values_to = "n") |>
  mutate(pct = 100 * n / n_peaks,
         motif = factor(motif, levels = c("GBS motif", "GAS Per box")))
message("events with motif: ", paste(sprintf("%s %s %d/%d", cmp$strain, cmp$motif,
                                             cmp$n, cmp$n_peaks), collapse = " | "))

p_cmp <- ggplot(cmp, aes(strain, pct, fill = motif)) +
  geom_col(position = position_dodge(width = .78), width = .66) +
  geom_text(aes(label = sprintf("%d/%d", n, n_peaks)),
            position = position_dodge(width = .78), vjust = -0.4,
            size = BAR_LABEL, fontface = "bold", colour = INK) +
  scale_fill_manual(values = MCOL) +
  scale_y_continuous(expand = expansion(mult = c(0, .2)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Binding events with motif") +
  theme_pub()
save_fig(p_cmp, "motif_gbs_vs_gas")

## ---------------------------------------------------------------------------
## 4. genome-wide occurrences: bound vs not bound
## ---------------------------------------------------------------------------
occ <- bind_rows(lapply(STRAINS, function(s) {
  gw <- suppressMessages(read_tsv(fimo_gbs[[s]], comment = "#", show_col_types = FALSE))
  gw <- gw[!is.na(gw$start), ]
  d  <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  bound <- sum(d$gbs_perbox_in_window == TRUE, na.rm = TRUE)
  tibble(strain = s, state = c("Bound", "Not bound"),
         n = c(bound, nrow(gw) - bound))
})) |>
  mutate(state = factor(state, levels = c("Not bound", "Bound")))
tot <- occ |> group_by(strain) |> summarise(t = sum(n), .groups = "drop")
message("occupancy: ", paste(sprintf("%s %d/%d bound",
        occ$strain[occ$state == "Bound"], occ$n[occ$state == "Bound"], tot$t),
        collapse = " | "))

p_occ <- ggplot(occ, aes(strain, n, fill = state)) +
  geom_col(width = .62) +
  geom_text(data = tot, aes(strain, t, label = t), inherit.aes = FALSE,
            vjust = -0.4, size = BAR_LABEL, fontface = "bold", colour = INK) +
  scale_fill_manual(values = c("Not bound" = GREY, "Bound" = TEAL)) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = "Genome-wide motif occurrences") +
  theme_pub()
save_fig(p_occ, "motif_occupancy")

sink(type = "message"); close(lg)
