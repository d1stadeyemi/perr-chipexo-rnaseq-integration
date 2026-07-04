## ============================================================================
## run_chipexoqual.R — ChIP-exo-specific QC for one BAM
## ============================================================================
## ChIPexoQual (Welch et al., NAR 2017) summarises three ChIP-exo QC axes:
##   * ARC vs URC plot  -> enrichment vs library-complexity structure
##   * FSR distribution -> forward-strand ratio (strand balance)
##   * beta1 / beta2    -> regression coefficients: beta1 tracks enrichment,
##                         beta2 tracks (inverse) library complexity; large
##                         beta2 flags low-complexity libraries.
## Minimum input is the aligned BAM (ExoData accepts a bam path directly).
## ============================================================================

suppressMessages({
  library(ChIPexoQual)
  library(GenomicAlignments)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

bam        <- snakemake@input[["bam"]]
sample_id  <- snakemake@params[["sample_id"]]
kind       <- snakemake@params[["kind"]]
strain     <- snakemake@params[["strain"]]
mc_cores   <- as.integer(snakemake@params[["mc_cores"]])

out_arc    <- snakemake@output[["arc_urc"]]
out_fsr    <- snakemake@output[["fsr"]]
out_rc     <- snakemake@output[["regioncomp"]]
out_pd     <- snakemake@output[["paramdist"]]
out_stats  <- snakemake@output[["stats"]]

message("ChIPexoQual on: ", bam)

## Build the ExoData object (islands + per-island stats) directly from the BAM.
ex <- ExoData(file = bam, mc.cores = mc_cores, verbose = TRUE)

## ---- diagnostic plots ------------------------------------------------------
save_plot <- function(p, path, w = 6, h = 5) {
  pdf(path, width = w, height = h)
  print(p)
  dev.off()
}

save_plot(ARCvURCplot(ex, names.input = sample_id), out_arc)
save_plot(FSRDistplot(ex, names.input = sample_id), out_fsr)
save_plot(regionCompplot(ex, names.input = sample_id), out_rc)
save_plot(paramDistBoxplot(ex, names.input = sample_id), out_pd)

## ---- summary metrics -------------------------------------------------------
## beta1 / beta2 are distributions (from repeated island subsampling); we
## report their medians. nreads = total reads used.
b1 <- tryCatch(beta1(ex), error = function(e) NA_real_)
b2 <- tryCatch(beta2(ex), error = function(e) NA_real_)
nr <- tryCatch(nreads(ex), error = function(e) NA_real_)

## Genome-wide forward-strand ratio (median over islands) for a strand-balance
## number that complements the FSR distribution plot.
fsr_median <- tryCatch(
  median(S4Vectors::mcols(ex)$FSR, na.rm = TRUE),
  error = function(e) NA_real_
)

med <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

stats <- data.frame(
  sample       = sample_id,
  strain       = strain,
  kind         = kind,
  n_reads      = nr,
  beta1_median = med(b1),   # enrichment (higher = better)
  beta2_median = med(b2),   # low-complexity flag (higher = worse)
  FSR_median   = fsr_median,
  n_islands    = length(ex),
  stringsAsFactors = FALSE
)

write.table(stats, file = out_stats, sep = "\t",
            quote = FALSE, row.names = FALSE)

message("wrote metrics: ", out_stats)
sink(type = "message"); sink(type = "output"); close(log)
