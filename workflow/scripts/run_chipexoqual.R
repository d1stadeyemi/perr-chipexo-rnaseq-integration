## ============================================================================
## run_chipexoqual.R — ChIP-exo QC for one BAM (small-genome safe, v2)
## ============================================================================
## Root cause of the earlier empty CJB111 panels:
##   ChIPexoQual's ExoData() samples `nregions` islands WITHOUT replacement to
##   fit the depth ~ uniquePos + width regression (beta1/beta2). On a ~2 Mb
##   genome there can be fewer islands than the default nregions=1000, which
##   throws "cannot take a sample larger than the population".
##   `nregions` ONLY controls that regression subsample; island formation (and
##   thus the ARC/URC, FSR, region-composition island stats) is independent.
##
## Fix vs v1:
##   v1 built a capped second object and PLOTTED THAT. When the capped build
##   degraded (CJB111), all three diagnostic plots emptied and beta2 -> None.
##   v2 plots the three island-level diagnostics from `ex_full`, built with
##   nregions=1 (can NEVER crash, still forms ALL islands). beta1/beta2 and the
##   paramDist boxplot come from a SEPARATE best-effort capped build; if that
##   fails, only the beta numbers degrade to NA — the diagnostic panels are
##   always complete.
## ============================================================================

suppressMessages({
  library(ChIPexoQual)
  library(GenomicAlignments)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

bam       <- snakemake@input[["bam"]]
sample_id <- snakemake@params[["sample_id"]]
kind      <- snakemake@params[["kind"]]
strain    <- snakemake@params[["strain"]]
mc_cores  <- as.integer(snakemake@params[["mc_cores"]])

out_arc   <- snakemake@output[["arc_urc"]]
out_fsr   <- snakemake@output[["fsr"]]
out_rc    <- snakemake@output[["regioncomp"]]
out_pd    <- snakemake@output[["paramdist"]]
out_stats <- snakemake@output[["stats"]]

message("ChIPexoQual on: ", bam)

save_plot <- function(expr, path, w = 6, h = 5) {
  pdf(path, width = w, height = h)
  tryCatch(print(expr), error = function(e) {
    plot.new(); title(main = paste("plot unavailable:", conditionMessage(e)))
  })
  dev.off()
}

## ---- ex_full: guaranteed-complete object for the island-level diagnostics --
## nregions=1 samples exactly 1 island for the (unused-here) regression, so it
## can never exceed the population; ALL islands are still formed.
ex_full <- ExoData(file = bam, mc.cores = mc_cores,
                   nregions = 1L, ntimes = 1L, verbose = TRUE)
n_islands <- length(ex_full)
message("n_islands = ", n_islands)

save_plot(ARCvURCplot(ex_full,   names.input = sample_id), out_arc)
save_plot(FSRDistplot(ex_full,   names.input = sample_id), out_fsr)
save_plot(regionCompplot(ex_full, names.input = sample_id), out_rc)

## ---- ex_param: best-effort regression for beta1/beta2 + paramDist boxplot --
nreg <- max(1L, min(1000L, n_islands))
ntim <- 100L
ex_param <- tryCatch(
  ExoData(file = bam, mc.cores = mc_cores,
          nregions = nreg, ntimes = ntim, verbose = FALSE),
  error = function(e) { message("param pass failed: ", conditionMessage(e)); NULL }
)

b1 <- b2 <- NA_real_
if (!is.null(ex_param)) {
  save_plot(paramDistBoxplot(ex_param, names.input = sample_id), out_pd)
  b1 <- tryCatch(median(beta1(ex_param), na.rm = TRUE), error = function(e) NA_real_)
  b2 <- tryCatch(median(beta2(ex_param), na.rm = TRUE), error = function(e) NA_real_)
} else {
  pdf(out_pd); plot.new()
  title(main = "paramDist unavailable (regression pass failed)"); dev.off()
}

fsr_median <- tryCatch(median(S4Vectors::mcols(ex_full)$FSR, na.rm = TRUE),
                       error = function(e) NA_real_)
nr <- tryCatch(nreads(ex_full), error = function(e)
        tryCatch(sum(S4Vectors::mcols(ex_full)$depth), error = function(e2) NA_real_))

stats <- data.frame(
  sample = sample_id, strain = strain, kind = kind,
  n_reads = nr, beta1_median = b1, beta2_median = b2,
  FSR_median = fsr_median, n_islands = n_islands,
  stringsAsFactors = FALSE
)
write.table(stats, file = out_stats, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote metrics: ", out_stats)
sink(type = "message"); sink(type = "output"); close(log)
