## ============================================================================
## plot_orthology.R — cross-serotype PerR regulon figures (ggplot2 + patchwork)
## ============================================================================
## Layout: 1 row x 2 panels per figure.
##
## READING CAUTION baked into the subtitles: CJB111 yielded ~4x fewer binding
## events, partly because it was prepared with half the standard bead/antibody
## slurry at lower depth. Strain-specific categories are therefore confounded
## with SENSITIVITY. SHARED_DIRECT is the robust category - both strains had to
## detect binding AND differential expression independently.
## ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(ggrepel); library(forcats)
})

xs_file <- snakemake@input[["crossstrain"]]
sm_file <- snakemake@input[["summary"]]
sa      <- snakemake@params[["strain_a"]]
sb      <- snakemake@params[["strain_b"]]
outdir  <- snakemake@params[["outdir"]]
logp    <- snakemake@log[[1]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt")

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"
GREY <- "#B4BEC7"; MUTE <- "#6B7A8A"; SAND <- "#C9743F"

theme_perr <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#EDF1F3", linewidth = .3),
          axis.text  = element_text(colour = INK),
          axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", hjust = 0),
          plot.subtitle = element_text(colour = MUTE, size = 8.5, hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(8, 12, 8, 8))
}

xs <- suppressMessages(read_tsv(xs_file, show_col_types = FALSE))
names(xs) <- gsub("^A909", sa, names(xs))   # tolerate literal names
ca <- paste0(sa, "_class");  cb <- paste0(sb, "_class")
fa <- paste0(sa, "_log2FC"); fb <- paste0(sb, "_log2FC")

## friendlier labels; OTHER is really "DE in one strain only, bound in neither"
relab <- c(SHARED_DIRECT       = "Shared direct",
           SHARED_BOUND        = "Shared bound",
           SHARED_INDIRECT     = "Shared indirect",
           DIVERGENT_BINDING   = "Divergent binding",
           OTHER               = "Unclassified")
relab[paste0(sa, "_ONLY_DIRECT")] <- paste0(sa, " only direct")
relab[paste0(sb, "_ONLY_DIRECT")] <- paste0(sb, " only direct")
relab[paste0("DE_", sa, "_ONLY")] <- paste0("DE ", sa, " only")
relab[paste0("DE_", sb, "_ONLY")] <- paste0("DE ", sb, " only")
relab[paste0("NO_ORTHOLOG_", sa)] <- paste0("No orthologue (", sa, ")")
relab[paste0("NO_ORTHOLOG_", sb)] <- paste0("No orthologue (", sb, ")")

PAL <- c("Shared direct" = RUST, "Shared bound" = TEAL,
         "Shared indirect" = "#7FB3B1", "Divergent binding" = SAND,
         "Unclassified" = "#DCE2E6")
PAL[paste0(sa, " only direct")]      <- "#D08C5A"
PAL[paste0(sb, " only direct")]      <- "#8C6BB1"
PAL[paste0("DE ", sa, " only")]      <- GREY
PAL[paste0("DE ", sb, " only")]      <- "#9AA7B1"
PAL[paste0("No orthologue (", sa, ")")] <- "#5C6B7A"
PAL[paste0("No orthologue (", sb, ")")] <- "#7D8B99"

xs <- xs %>%
  mutate(cons_lab = ifelse(conservation %in% names(relab),
                           relab[conservation], conservation))

## ---------------------------------------------------------------- FIGURE 1
## P1: conservation composition
cnt <- xs %>% count(cons_lab, name = "n") %>% arrange(n) %>%
  mutate(cons_lab = fct_inorder(cons_lab),
         robust = cons_lab == "Shared direct")

p1 <- ggplot(cnt, aes(cons_lab, n, fill = cons_lab)) +
  geom_col(width = .68) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.1, colour = INK) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "Conservation of the PerR regulon",
       subtitle = paste0("1:1 orthologue pairs (RBH), plus genes with no ",
                         "orthologue. Strain-specific classes are confounded ",
                         "with sensitivity."),
       x = NULL, y = "orthologue pairs") +
  theme_perr()

## P2: cross-strain expression scatter
sc <- xs %>% filter(!is.na(.data[[fa]]), !is.na(.data[[fb]]))
lab_sc <- sc %>% filter(conservation == "SHARED_DIRECT") %>%
  mutate(lab = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name,
                      .data[[paste0(sa, "_gene")]]))
r <- if (nrow(sc) > 2) cor(sc[[fa]], sc[[fb]], use = "complete.obs") else NA_real_

p2 <- ggplot(sc, aes(.data[[fa]], .data[[fb]])) +
  geom_hline(yintercept = 0, colour = "#D8DEE3", linewidth = .4) +
  geom_vline(xintercept = 0, colour = "#D8DEE3", linewidth = .4) +
  geom_abline(slope = 1, intercept = 0, colour = "#D8DEE3",
              linetype = "dashed", linewidth = .4) +
  geom_point(aes(colour = cons_lab, size = cons_lab == "Shared direct"),
             alpha = .8) +
  ggrepel::geom_text_repel(data = lab_sc, aes(label = lab), size = 2.7,
                           colour = RUST, max.overlaps = 25,
                           segment.colour = "#D8DEE3") +
  scale_colour_manual(values = PAL) +
  scale_size_manual(values = c(`FALSE` = 1.3, `TRUE` = 2.8), guide = "none") +
  labs(title = "Expression response across serotypes",
       subtitle = sprintf("Pearson r = %.2f; dashed line = identity", r),
       x = paste0(sa, "  log2 fold change"),
       y = paste0(sb, "  log2 fold change")) +
  theme_perr() +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 2.4)))

## ---------------------------------------------------------------- FIGURE 0
## P0: binding conservation (expression-independent).
## This is the axis the half-slurry caveat acts on. The paired panel P1 mixes
## binding with expression response and therefore UNDERSTATES conservation of
## the cistrome; both are shown so the distinction is explicit.
bind_lab <- c(BOUND_BOTH = "Bound in both", BOUND_NEITHER = "Bound in neither")
bind_lab[paste0("BOUND_", sa)] <- paste0("Bound ", sa, " only")
bind_lab[paste0("BOUND_", sb)] <- paste0("Bound ", sb, " only")

BPAL <- c("Bound in both" = RUST, "Bound in neither" = "#DCE2E6")
BPAL[paste0("Bound ", sa, " only")] <- "#D08C5A"
BPAL[paste0("Bound ", sb, " only")] <- "#8C6BB1"

## Restrict to genes present in BOTH proteomes. A gene with no orthologue
## cannot be "bound in both", so including it would inflate the strain-specific
## bars with a gene-CONTENT difference rather than a BINDING difference.
bcnt <- xs %>%
  filter(binding != "BOUND_NEITHER",
         nzchar(.data[[paste0(sa, "_gene")]]),
         nzchar(.data[[paste0(sb, "_gene")]])) %>%
  mutate(bind_lab = ifelse(binding %in% names(bind_lab),
                           bind_lab[binding], binding)) %>%
  count(bind_lab, name = "n") %>% arrange(n) %>%
  mutate(bind_lab = fct_inorder(bind_lab))

p0 <- ggplot(bcnt, aes(bind_lab, n, fill = bind_lab)) +
  geom_col(width = .68) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.1, colour = INK) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = BPAL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(title = "Conservation of PerR BINDING",
       subtitle = paste0("binding only -- no expression threshold applied; ",
                         "strain-specific binding is confounded with sensitivity"),
       x = NULL, y = "orthologue pairs") +
  theme_perr()

fig0 <- (p0 | p1) +
  plot_annotation(
    title = "Binding conservation vs binding x expression conservation",
    subtitle = paste0("left: cistrome only.  right: binding combined with the ",
                      "RNA-seq response.  The right panel is the stricter test."),
    theme = theme(plot.title = element_text(face = "bold", size = 14, colour = INK),
                  plot.subtitle = element_text(colour = MUTE, size = 9)))
ggsave(file.path(outdir, "orthology_fig0_binding_axis.pdf"), fig0,
       width = 13, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, "orthology_fig0_binding_axis.png"), fig0,
       width = 13, height = 5.6, dpi = 300)

fig1 <- (p1 | p2) +
  plot_annotation(title = "Cross-serotype comparison of the PerR regulon",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14, colour = INK)))
ggsave(file.path(outdir, "orthology_fig1_conservation.pdf"), fig1,
       width = 13, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, "orthology_fig1_conservation.png"), fig1,
       width = 13, height = 5.6, dpi = 300)

## ---------------------------------------------------------------- FIGURE 2
## P3: the conserved core - dumbbell of both strains' fold changes
core <- xs %>% filter(conservation == "SHARED_DIRECT") %>%
  mutate(lab = ifelse(!is.na(gene_name) & nzchar(gene_name),
                      paste0(.data[[paste0(sa, "_gene")]], " (", gene_name, ")"),
                      .data[[paste0(sa, "_gene")]])) %>%
  arrange(.data[[fa]]) %>% mutate(lab = fct_inorder(lab))

if (nrow(core)) {
  core_long <- core %>%
    select(lab, !!sa := all_of(fa), !!sb := all_of(fb)) %>%
    pivot_longer(-lab, names_to = "strain", values_to = "log2FC")
  p3 <- ggplot(core, aes(y = lab)) +
    geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
    geom_segment(aes(x = .data[[fa]], xend = .data[[fb]], yend = lab),
                 colour = "#C7D0D6", linewidth = 1.1) +
    geom_point(data = core_long, aes(x = log2FC, colour = strain), size = 3.2) +
    scale_colour_manual(values = setNames(c(RUST, TEAL), c(sa, sb))) +
    labs(title = "Conserved core: direct in both serotypes",
         subtitle = "each line joins the two serotypes' expression response",
         x = "RNA-seq log2 fold change", y = NULL) +
    theme_perr()
} else {
  p3 <- ggplot() + theme_void() +
    annotate("text", 0, 0, label = "no shared direct targets", colour = MUTE)
}

## P4: class concordance matrix
mat <- xs %>%
  filter(!is.na(.data[[ca]]), !is.na(.data[[cb]]),
         nzchar(.data[[ca]]), nzchar(.data[[cb]])) %>%
  count(.data[[ca]], .data[[cb]], name = "n")
names(mat)[1:2] <- c("A", "B")

p4 <- ggplot(mat, aes(A, B, fill = n)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = n), size = 3, colour = INK) +
  scale_fill_gradient(low = "#EAF2F1", high = TEAL, guide = "none") +
  labs(title = "Classification concordance",
       subtitle = "orthologue pairs classified in both serotypes",
       x = paste0(sa, " classification"), y = paste0(sb, " classification")) +
  theme_perr() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig2 <- (p3 | p4) +
  plot_annotation(title = "Conserved core and classification concordance",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14, colour = INK)))
ggsave(file.path(outdir, "orthology_fig2_core_concordance.pdf"), fig2,
       width = 13, height = 5.8, device = cairo_pdf)
ggsave(file.path(outdir, "orthology_fig2_core_concordance.png"), fig2,
       width = 13, height = 5.8, dpi = 300)

writeLines(sprintf("rows total: %d", nrow(xs)), lg)
writeLines(sprintf("1:1 orthologue pairs: %d",
                   sum(!grepl("^NO_ORTHOLOG", xs$conservation))), lg)
writeLines(sprintf("shared direct: %d", nrow(core)), lg)
writeLines(sprintf("cross-strain log2FC correlation: %.3f", r), lg)
writeLines(sprintf("bound in both strains: %d",
                   sum(xs$binding == "BOUND_BOTH")), lg)
writeLines(sprintf("no-orthologue rows: %d",
                   sum(grepl("^NO_ORTHOLOG", xs$conservation))), lg)
close(lg)