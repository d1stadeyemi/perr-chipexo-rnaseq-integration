## =============================================================================
## plot_motif_evidence.R -- the sequence evidence for the GBS PerR motif
## =============================================================================
## Four panels answering four separate questions, none of which depends on a
## per-site FIMO q-value (those run 0.6-1.0 for a genome-wide 11-mer scan and
## must never be quoted):
##
##   A  what does the motif look like            -> sequence logo from the PWM
##   B  is it AT the crosslink or merely nearby   -> offset distribution
##   C  does it beat the S. pyogenes consensus    -> like-for-like peak counts
##   D  does motif presence predict binding       -> genome-wide vs under peaks
##
## Panel B is the strongest evidence in the study: a median offset of ~0 bp
## cannot arise from a scanning artefact.
## Panel D is the argument for why ChIP-exo was necessary rather than a scan.
## =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(ggseqlogo); library(scales)
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
on.exit({sink(type = "message"); close(lg)}, add = TRUE)

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"
MUTE <- "#6B7A8A"; GREY <- "#B4BEC7"; PALE <- "#EDF1F3"
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
## PANEL A -- sequence logo, read straight from the MEME PWM that FIMO scanned
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
consensus <- paste(c("A","C","G","T")[apply(pwm, 2, which.max)], collapse = "")
message("PWM width ", ncol(pwm), "  majority consensus ", consensus)

pA <- ggseqlogo(pwm, method = "bits") +
  labs(title = "A   The GBS PerR motif",
       subtitle = paste0("de novo, MEME; ", ncol(pwm),
                         " bp; one matrix scans both genomes"),
       x = "position", y = "bits") +
  theme_perr() + theme(legend.position = "none")

## ---------------------------------------------------------------------------
## PANEL B -- offset from the crosslink
## Every hit in the 500 bp discovery window is plotted, not only those that
## passed the distance filter, so the filter's effect is visible rather than
## assumed. A scanning artefact would be flat across the window.
## ---------------------------------------------------------------------------
off <- bind_rows(lapply(STRAINS, function(s) {
  d <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  tibble(strain = s,
         offset = suppressWarnings(as.numeric(d$gbs_best_offset))) |>
    filter(!is.na(offset))
}))
med <- off |> group_by(strain) |> summarise(m = median(offset), n = n(), .groups = "drop")
for (k in seq_len(nrow(med)))
  message(sprintf("%s: %d motif hits in window, median offset %+d bp",
                  med$strain[k], med$n[k], as.integer(med$m[k])))

pB <- ggplot(off, aes(offset, fill = strain)) +
  annotate("rect", xmin = -HALFWIN, xmax = HALFWIN, ymin = -Inf, ymax = Inf,
           fill = PALE, alpha = .9) +
  geom_histogram(binwidth = 10, colour = "white", linewidth = .2, alpha = .85) +
  geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
  geom_vline(data = med, aes(xintercept = m, colour = strain),
             linetype = "dashed", linewidth = .6, show.legend = FALSE) +
  facet_wrap(~ strain, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = SCOL, guide = "none") +
  scale_colour_manual(values = SCOL, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks(6)) +
  labs(title = "B   The motif sits ON the crosslink",
       subtitle = paste0("shaded = the +/-", HALFWIN,
                         " bp window used to call a peak motif-positive;",
                         "  dashed = median"),
       x = "motif centre relative to the ChExMix crosslink (bp)",
       y = "motif hits") +
  theme_perr()

## ---------------------------------------------------------------------------
## PANEL C -- GAS box vs GBS motif, judged by the identical rule
## Both scans use the same threshold, the same background and the same distance
## cutoff, so the comparison is like for like.
## ---------------------------------------------------------------------------
cmp <- bind_rows(lapply(STRAINS, function(s) {
  d <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  tibble(strain = s, n_peaks = nrow(d),
         GBS = sum(d$gbs_perbox_in_window == TRUE, na.rm = TRUE),
         GAS = sum(d$perbox_in_window     == TRUE, na.rm = TRUE))
})) |>
  pivot_longer(c(GBS, GAS), names_to = "motif", values_to = "n") |>
  mutate(pct = 100 * n / n_peaks,
         motif = factor(motif, levels = c("GAS", "GBS"),
                        labels = c("S. pyogenes Per box", "GBS de novo motif")))
message("peaks with motif: ", paste(sprintf("%s %s %d/%d", cmp$strain, cmp$motif,
                                            cmp$n, cmp$n_peaks), collapse = " | "))

pC <- ggplot(cmp, aes(strain, pct, fill = motif)) +
  geom_col(position = position_dodge(width = .75), width = .62) +
  geom_text(aes(label = sprintf("%d/%d", n, n_peaks)),
            position = position_dodge(width = .75), vjust = -0.4,
            size = 3, colour = INK) +
  scale_fill_manual(values = c("S. pyogenes Per box" = GREY,
                               "GBS de novo motif" = RUST)) +
  scale_y_continuous(expand = expansion(mult = c(0, .18)),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "C   The GBS motif explains binding; the GAS box does not",
       subtitle = "same FIMO threshold, same background, same distance rule",
       x = NULL, y = "peaks containing the motif") +
  theme_perr()

## ---------------------------------------------------------------------------
## PANEL D -- occupancy
## Most genome-wide motif occurrences are NOT bound. This is the reason binding
## data was needed rather than a consensus scan.
## ---------------------------------------------------------------------------
occ <- bind_rows(lapply(STRAINS, function(s) {
  gw <- suppressMessages(read_tsv(fimo_gbs[[s]], comment = "#", show_col_types = FALSE))
  gw <- gw[!is.na(gw$start), ]
  d  <- suppressMessages(read_tsv(class_files[[s]], show_col_types = FALSE))
  bound <- sum(d$gbs_perbox_in_window == TRUE, na.rm = TRUE)
  tibble(strain = s,
         state = c("under a peak", "not bound"),
         n = c(bound, nrow(gw) - bound))
})) |>
  mutate(state = factor(state, levels = c("not bound", "under a peak")))
tot <- occ |> group_by(strain) |> summarise(t = sum(n), .groups = "drop")
message("occupancy: ", paste(sprintf("%s %d/%d bound",
        occ$strain[occ$state == "under a peak"],
        occ$n[occ$state == "under a peak"], tot$t), collapse = " | "))

pD <- ggplot(occ, aes(strain, n, fill = state)) +
  geom_col(width = .62) +
  geom_text(data = occ |> filter(state == "under a peak"),
            aes(label = sprintf("%d bound", n)), vjust = 1.6,
            size = 3, colour = "white") +
  geom_text(data = tot, aes(strain, t, label = sprintf("%d motif sites", t)),
            inherit.aes = FALSE, vjust = -0.5, size = 3, colour = INK) +
  scale_fill_manual(values = c("not bound" = GREY, "under a peak" = TEAL)) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "D   Most motif occurrences are not bound",
       subtitle = "a consensus scan predicts far more sites than PerR occupies",
       x = NULL, y = "genome-wide motif occurrences") +
  theme_perr()

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Sequence evidence for the GBS PerR binding motif",
    subtitle = "Discovered de novo in both serotypes; distinct from the S. pyogenes Per box",
    theme = theme(plot.title = element_text(face = "bold", size = 15, colour = INK),
                  plot.subtitle = element_text(colour = MUTE, size = 10)))

ggsave(file.path(outdir, "motif_evidence.pdf"), fig, width = 12, height = 9)
ggsave(file.path(outdir, "motif_evidence.png"), fig, width = 12, height = 9, dpi = 300)
message("wrote motif_evidence.pdf")