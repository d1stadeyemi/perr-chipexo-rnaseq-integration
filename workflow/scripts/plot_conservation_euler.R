## =============================================================================
## plot_conservation_euler.R -- cross-serotype binding and regulon overlap
## =============================================================================
## Area-proportional Euler rather than a fixed-circle Venn, so the areas mean
## what they appear to mean.
##
## IMPORTANT: the A909-only region is NOT a clean biological result. The
## sequencing core prepared the CJB111 samples with half the standard bead and
## antibody slurry, and the orthologous-position test (module 09) estimates that
## a majority of A909-only sites carry sub-threshold signal at the corresponding
## CJB111 position. The caption states this on the figure itself, because a Venn
## is read at a glance and would otherwise be remembered as a 20-fold biological
## difference.
## =============================================================================
suppressPackageStartupMessages({
  library(eulerr); library(dplyr); library(readr); library(grid)
})

cross   <- snakemake@input[["cross"]]
summ    <- snakemake@input[["summary"]]
outdir  <- snakemake@params[["outdir"]]
strains <- unlist(snakemake@params[["strains"]])
logp    <- snakemake@log[[1]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt"); sink(lg, type = "message")
on.exit({sink(type = "message"); close(lg)}, add = TRUE)

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"; MUTE <- "#6B7A8A"
A <- strains[1]; B <- strains[2]

## Read the categories the pipeline already assigned rather than recomputing
## them. crossstrain_regulon.py applies the orthologue and biotype rules that a
## naive startsWith() on the classification columns would miss, and the two
## must not be allowed to disagree.
s <- suppressMessages(read_tsv(summ, show_col_types = FALSE))
gv <- function(cat) {
  v <- s$n[s$category == cat]
  if (!length(v)) stop("category '", cat, "' not in ", summ)
  as.integer(v[1])
}
bind <- c(setNames(gv(paste0("BOUND_", A)), A),
          setNames(gv(paste0("BOUND_", B)), B),
          setNames(gv("BOUND_BOTH"), paste(A, B, sep = "&")))
dirc <- c(setNames(gv(paste0(A, "_ONLY_DIRECT")), A),
          setNames(gv(paste0(B, "_ONLY_DIRECT")), B),
          setNames(gv("SHARED_DIRECT"), paste(A, B, sep = "&")))
message("binding : ", paste(names(bind), bind, collapse = "  "))
message("regulon : ", paste(names(dirc), dirc, collapse = "  "))

draw <- function(v, ttl, sub) {
  plot(euler(v, shape = "ellipse"),
       fills = list(fill = c(RUST, TEAL), alpha = .45),
       edges = list(col = c(RUST, TEAL), lwd = 2),
       labels = list(col = INK, fontfamily = "sans", cex = 1.0),
       quantities = list(col = INK, fontfamily = "sans", cex = 1.1),
       main = list(label = ttl, col = INK, fontface = "bold", cex = 1.15))
}

pdf(file.path(outdir, "conservation_euler.pdf"), width = 11, height = 5.6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(1, 3.2), "null"))))
grid.text("Cross-serotype conservation of PerR binding",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .72, gp = gpar(fontface = "bold", cex = 1.4, col = INK))
grid.text(paste0("left: bound in either serotype (no expression threshold).   ",
                 "right: bound AND differentially expressed."),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .34, gp = gpar(cex = .85, col = MUTE))
grid.text(paste0("The ", A, "-only region is confounded: CJB111 was prepared with half the ",
                 "standard bead and antibody slurry, and the orthologous-position test finds ",
                 "sub-threshold signal at most of these sites."),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .06, gp = gpar(cex = .75, col = MUTE, fontface = "italic"))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
print(draw(bind, "PerR binding (cistrome)", ""), newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
print(draw(dirc, "Direct regulon (binding + expression)", ""), newpage = FALSE); popViewport()
dev.off()

png(file.path(outdir, "conservation_euler.png"), width = 11, height = 5.6,
    units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(1, 3.2), "null"))))
grid.text("Cross-serotype conservation of PerR binding",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .72, gp = gpar(fontface = "bold", cex = 1.4, col = INK))
grid.text(paste0("left: bound in either serotype (no expression threshold).   ",
                 "right: bound AND differentially expressed."),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .34, gp = gpar(cex = .85, col = MUTE))
grid.text(paste0("The ", A, "-only region is confounded: CJB111 was prepared with half the ",
                 "standard bead and antibody slurry, and the orthologous-position test finds ",
                 "sub-threshold signal at most of these sites."),
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
          y = .06, gp = gpar(cex = .75, col = MUTE, fontface = "italic"))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
print(draw(bind, "PerR binding (cistrome)", ""), newpage = FALSE); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
print(draw(dirc, "Direct regulon (binding + expression)", ""), newpage = FALSE); popViewport()
dev.off()
message("wrote conservation_euler.pdf")