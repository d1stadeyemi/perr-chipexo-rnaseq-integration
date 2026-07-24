## ============================================================================
## plot_regulon.R — ggplot2 figures for the ChIP-exo x RNA-seq direct regulon
## ============================================================================
## Layout: each figure is 1 row x 2 panels (patchwork), so four plots become
## two figures per strain, plus a standalone cross-strain composition bar.
##
## COUNTING: the regulon table has one row per (peak, gene) pair, so a gene
## bound by several peaks appears several times. Everything here is collapsed
## to UNIQUE GENES using a class precedence
##   DIRECT_PERBOX > DIRECT > BOUND_NOT_DE > BOUND_NO_RNA > INDIRECT
## so a gene that is DIRECT at one peak is never also counted elsewhere.
##
## BOUND_NO_RNA is a DATA-COMPLETENESS artefact (a bound gene whose locus tag
## has no row in the RNA-seq table), not a biological category. It is retained
## in the log and in regulon_gene_counts.tsv for transparency but EXCLUDED from
## every figure.
## ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(ggrepel)
})

reg_files <- unlist(snakemake@input[["regulon"]])
rna_files <- unlist(snakemake@input[["rnaseq"]])
strains   <- unlist(snakemake@params[["strains"]])
outdir    <- snakemake@params[["outdir"]]
logp      <- snakemake@log[[1]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt")

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"; GREY <- "#B4BEC7"
PREC <- c("DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "BOUND_NO_RNA", "INDIRECT")
## classes shown in figures (BOUND_NO_RNA is an artefact - see header)
PLOT_CLASSES <- c("DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "INDIRECT")
PAL  <- c(DIRECT_PERBOX = RUST, DIRECT = "#C9743F", BOUND_NOT_DE = TEAL,
          BOUND_NO_RNA = "#7FB3B1", INDIRECT = GREY)

theme_perr <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#EDF1F3", linewidth = .3),
          axis.text  = element_text(colour = INK),
          axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold", hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(8, 12, 8, 8))
}

## collapse peak-gene rows to one row per gene (best class wins)
gene_level <- function(reg) {
  reg %>%
    mutate(.rank = match(classification, PREC)) %>%
    arrange(.rank) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(-.rank)
}

all_counts <- list()

for (i in seq_along(strains)) {
  strain <- strains[i]
  reg <- suppressMessages(read_tsv(reg_files[i], col_types = cols(.default = "c"))) %>%
    mutate(across(c(rna_log2FC, rna_padj, chip_log2FC, rna_baseMean),
                  ~suppressWarnings(as.numeric(.))))
  rna <- suppressMessages(read_tsv(rna_files[i], show_col_types = FALSE))

  g <- gene_level(reg)
  cnt_all <- g %>% count(classification, name = "n_genes")
  all_counts[[strain]] <- cnt_all %>% mutate(strain = strain)   # full, for the TSV
  cnt <- cnt_all %>%
    filter(classification %in% PLOT_CLASSES) %>%
    mutate(classification = factor(classification, levels = rev(PLOT_CLASSES))) %>%
    arrange(classification)

  writeLines(sprintf("%s: rows=%d unique_genes=%d", strain, nrow(reg), nrow(g)), lg)
  for (k in seq_len(nrow(cnt_all)))
    writeLines(sprintf("   %s: %d%s", cnt_all$classification[k], cnt_all$n_genes[k],
               ifelse(cnt_all$classification[k] == "BOUND_NO_RNA",
                      "   [artefact: excluded from figures]", "")), lg)

  direct <- g %>% filter(grepl("^DIRECT", classification))
  bound  <- g %>% filter(grepl("^(DIRECT|BOUND)", classification)) %>% pull(gene_id)

  ## ---- P1: composition -------------------------------------------------
  p1 <- ggplot(cnt, aes(classification, n_genes, fill = classification)) +
    geom_col(width = .68) +
    geom_text(aes(label = n_genes), hjust = -0.2, size = 3.2, colour = INK) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "Regulon composition", x = NULL, y = "unique genes") +
    theme_perr()

  ## ---- P2: volcano -----------------------------------------------------
  v <- rna %>%
    filter(!is.na(log2FC), !is.na(padj)) %>%
    mutate(nlp   = -log10(pmax(padj, 1e-300)),
           class = case_when(gene_id %in% direct$gene_id ~ "direct target",
                             gene_id %in% bound          ~ "PerR-bound",
                             TRUE                        ~ "not bound"),
           class = factor(class, levels = c("not bound", "PerR-bound", "direct target")))
  lab <- v %>% filter(class == "direct target") %>%
    left_join(direct %>% select(gene_id, gene_name), by = "gene_id") %>%
    mutate(lab = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_id))

  p2 <- ggplot(v, aes(log2FC, nlp, colour = class, size = class)) +
    geom_vline(xintercept = 0, colour = "#D8DEE3", linewidth = .4) +
    geom_point(data = ~subset(., class == "not bound"), alpha = .5) +
    geom_point(data = ~subset(., class != "not bound")) +
    ggrepel::geom_text_repel(data = lab, aes(label = lab), size = 2.6,
                             colour = RUST, max.overlaps = 30,
                             segment.colour = "#D8DEE3", show.legend = FALSE) +
    scale_colour_manual(values = c("not bound" = GREY, "PerR-bound" = TEAL,
                                   "direct target" = RUST)) +
    scale_size_manual(values = c("not bound" = .7, "PerR-bound" = 1.5,
                                 "direct target" = 2.4), guide = "none") +
    labs(title = "Expression landscape",
         x = "RNA-seq log2 fold change", y = expression(-log[10]~adjusted~p)) +
    theme_perr()

  fig1 <- (p1 | p2) +
    plot_annotation(title = paste0("PerR direct regulon - ", strain),
                    theme = theme(plot.title = element_text(face = "bold",
                                                            size = 14, colour = INK)))
  ggsave(file.path(outdir, sprintf("%s_fig1_composition_volcano.pdf", strain)),
         fig1, width = 12, height = 5.2)
  ggsave(file.path(outdir, sprintf("%s_fig1_composition_volcano.png", strain)),
         fig1, width = 12, height = 5.2, dpi = 300)

  ## ---- P3: binding strength vs expression ------------------------------
  bd <- reg %>% filter(relationship != "none",
                       !is.na(chip_log2FC), !is.na(rna_log2FC)) %>%
    mutate(grp = ifelse(grepl("^DIRECT", classification), "direct", "bound"))
  p3 <- ggplot(bd, aes(chip_log2FC, rna_log2FC, colour = grp)) +
    geom_hline(yintercept = 0, colour = "#D8DEE3", linewidth = .4) +
    geom_point(aes(size = grp), alpha = .75) +
    scale_colour_manual(values = c(bound = TEAL, direct = RUST)) +
    scale_size_manual(values = c(bound = 1.5, direct = 2.6), guide = "none") +
    labs(title = "Binding strength vs expression change",
         x = "ChIP-exo log2 fold (peak vs control)",
         y = "RNA-seq log2 fold change") +
    theme_perr()

  ## ---- P4: direct target profile ---------------------------------------
  dd <- direct %>%
    mutate(lab = ifelse(!is.na(gene_name) & nzchar(gene_name),
                        paste0(gene_id, " (", gene_name, ")"), gene_id)) %>%
    arrange(rna_log2FC) %>% mutate(lab = factor(lab, levels = lab))
  p4 <- ggplot(dd, aes(lab, rna_log2FC, fill = classification)) +
    geom_col(width = .7) +
    geom_hline(yintercept = 0, colour = INK, linewidth = .4) +
    coord_flip() +
    scale_fill_manual(values = PAL, labels = function(x) gsub("_", " ", x)) +
    labs(title = "Direct targets", x = NULL, y = "RNA-seq log2 fold change") +
    theme_perr()

  fig2 <- (p3 | p4) +
    plot_annotation(title = paste0("PerR direct regulon - ", strain),
                    theme = theme(plot.title = element_text(face = "bold",
                                                            size = 14, colour = INK)))
  ggsave(file.path(outdir, sprintf("%s_fig2_binding_targets.pdf", strain)),
         fig2, width = 12, height = 5.8)
  ggsave(file.path(outdir, sprintf("%s_fig2_binding_targets.png", strain)),
         fig2, width = 12, height = 5.8, dpi = 300)
  writeLines(sprintf("  wrote figures for %s", strain), lg)
}

## ---- cross-strain composition -------------------------------------------
cc_all <- bind_rows(all_counts) %>%
  mutate(classification = factor(classification, levels = PREC))
write_tsv(cc_all, file.path(outdir, "regulon_gene_counts.tsv"))   # full, incl. artefact

cc <- cc_all %>%
  filter(classification %in% PLOT_CLASSES) %>%
  mutate(classification = factor(classification, levels = PLOT_CLASSES))

## Grouped (dodged) bars: strain on the x-axis, one bar per class side by side.
## Value labels are printed because the categories differ by ~2 orders of
## magnitude, which would otherwise make the small classes unreadable.
pcc <- ggplot(cc, aes(strain, n_genes, fill = classification)) +
  geom_col(position = position_dodge(width = .8), width = .72) +
  geom_text(aes(label = n_genes),
            position = position_dodge(width = .8),
            vjust = -0.35, size = 3, colour = INK) +
  scale_fill_manual(values = PAL, labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "PerR regulon composition by strain",
       subtitle = "unique genes; grouped bars per class",
       x = NULL, y = "unique genes") +
  theme_perr() +
  theme(plot.subtitle = element_text(colour = "#6B7A8A", size = 9))
ggsave(file.path(outdir, "regulon_composition.pdf"), pcc, width = 8, height = 5)
ggsave(file.path(outdir, "regulon_composition.png"), pcc, width = 8, height = 5, dpi = 300)

writeLines("wrote cross-strain composition", lg)
close(lg)