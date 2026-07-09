## ============================================================================
## run_chipexoqual.R — ChIP-exo QC for one BAM (small-genome safe, v3)
## ============================================================================
## Why this looks the way it does:
##   On a ~2 Mb genome at full depth, every read overlaps its neighbour and the
##   whole chromosome merges into ONE island (measured: CJB111 full = 1 island,
##   5% = 311, 1% = 7593; A909 full = 648). ChIPexoQual cannot characterise a
##   1-island library, which is why the CJB111 panels came out empty.
##
##   Fix: reads are subsampled to a FIXED read count BEFORE this script runs
##   (rule subsample_bam_for_qc). A fixed count — not a fixed fraction — keeps
##   beta1/beta2 comparable across samples of differing depth.
##
##   Second trap: calculateParamDist() does `stats[sample(.N, nregions)]` with
##   no internal guard, so nregions MUST be <= n_islands. ExoData(nregions=1L)
##   never trips this and still forms ALL islands, so we use it for the object
##   and compute the regression ourselves with a capped nregions.
##   Model is lm(depth ~ 0 + uniquePos + width): beta1 = uniquePos coefficient,
##   beta2 = width coefficient.
## ============================================================================

suppressMessages({
  library(ChIPexoQual)
  library(GenomicAlignments)
  library(data.table)
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

message("ChIPexoQual on (subsampled BAM): ", bam)

save_plot <- function(expr, path, w = 6, h = 5) {
  pdf(path, width = w, height = h)
  tryCatch(print(expr), error = function(e) {
    plot.new(); title(main = paste("plot unavailable:", conditionMessage(e)))
  })
  dev.off()
}

## ---- build ExoData (nregions=1L is always safe; keeps all islands) ----------
ex <- ExoData(file = bam, mc.cores = mc_cores,
              nregions = 1L, ntimes = 1L, verbose = TRUE)
n_islands <- length(ex)
message("islands: ", n_islands)

## ---- island-level diagnostic panels ----------------------------------------
save_plot(ARCvURCplot(ex,    names.input = sample_id), out_arc)
save_plot(FSRDistplot(ex,    names.input = sample_id), out_fsr)
save_plot(regionCompplot(ex, names.input = sample_id), out_rc)

## ---- beta1 / beta2 with nregions strictly <= n_islands ---------------------
b1 <- b2 <- NA_real_
if (n_islands >= 2L) {
  DT   <- formatRegions(ex)
  nreg <- min(1000L, n_islands)
  message("regression nregions: ", nreg)
  pd <- tryCatch(
    rbindlist(lapply(seq_len(100L), calculateParamDist, stats = DT, nregions = nreg)),
    error = function(e) { message("paramDist failed: ", conditionMessage(e)); NULL })
  if (!is.null(pd)) {
    b1 <- median(pd[term == "uniquePos", estimate], na.rm = TRUE)
    b2 <- median(pd[term == "width",     estimate], na.rm = TRUE)
    pdf(out_pd, width = 6, height = 5)
    boxplot(pd[term == "uniquePos", estimate], pd[term == "width", estimate],
            names = c("beta1", "beta2"), main = sample_id,
            ylab = "coefficient estimate")
    dev.off()
  }
}
if (!file.exists(out_pd)) {
  pdf(out_pd); plot.new(); title("paramDist unavailable (too few islands)"); dev.off()
}

## ---- summary metrics --------------------------------------------------------
fsr_median <- tryCatch(median(S4Vectors::mcols(ex)$FSR, na.rm = TRUE),
                       error = function(e) NA_real_)
nr <- tryCatch(nreads(ex), error = function(e)
        tryCatch(sum(S4Vectors::mcols(ex)$depth), error = function(e2) NA_real_))

stats <- data.frame(
  sample = sample_id, strain = strain, kind = kind,
  n_reads_subsampled = nr,
  beta1_median = b1, beta2_median = b2,
  FSR_median = fsr_median, n_islands = n_islands,
  stringsAsFactors = FALSE
)
write.table(stats, file = out_stats, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote metrics: ", out_stats)
sink(type = "message"); sink(type = "output"); close(log)