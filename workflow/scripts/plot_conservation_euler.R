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
outdir  <- snakemake@params[["outdir"]]
strains <- unlist(snakemake@params[["strains"]])
logp    <- snakemake@log[[1]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(logp), recursive = TRUE, showWarnings = FALSE)
lg <- file(logp, open = "wt"); sink(lg, type = "message")
on.exit({sink(type = "message"); close(lg)}, add = TRUE)

INK <- "#1F2933"; TEAL <- "#0E6E6B"; RUST <- "#A8471C"; MUTE <- "#6B7A8A"
A <- strains[1]; B <- strains[2]

d <- suppressMessages(read_tsv(cross, show_col_types = FALSE))
is_bound  <- function(x) !is.na(x) & grepl("^(DIRECT|BOUND)", x)
is_direct <- function(x) !is.na(x) & grepl("^DIRECT", x)

ca <- d[[grep(paste0("^", A, ".*class"), names(d))[1]]]
cb <- d[[grep(paste0("^", B, ".*class"), names(d))[1]]]
if (is.null(ca) || is.null(cb)) {
  ca <- d[[5]]; cb <- d[[6]]
  message("NOTE: classification columns located positionally (5,6); check the header")
}

bind <- c(setNames(sum( is_bound(ca) & !is_bound(cb)), A),
          setNames(sum(!is_bound(ca) &  is_bound(cb)), B),
          setNames(sum( is_bound(ca) &  is_bound(cb)), paste(A, B, sep = "&")))
dirc <- c(setNames(sum( is_direct(ca) & !is_direct(cb)), A),
          setNames(sum(!is_direct(ca) &  is_direct(cb)), B),
          setNames(sum( is_direct(ca) &  is_direct(cb)), paste(A, B, sep = "&")))
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