## run_chipexoqual.R — ChIP-exo QC for one BAM (small-genome safe)
## Fix: ExoData() subsamples `nregions` islands WITHOUT replacement; on a ~2 Mb
## bacterial genome there can be fewer islands than the default 1000, which
## triggers "cannot take a sample larger than the population". We cap nregions
## at the number of islands actually present.
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

## Pass 1: build cheaply (nregions=1) just to count islands. This does NOT run
## the failing large subsample.
ex0 <- ExoData(file = bam, mc.cores = mc_cores,
               nregions = 1L, ntimes = 1L, verbose = FALSE)
n_islands <- length(ex0)
message("n_islands = ", n_islands)

## Pass 2: rebuild with a subsample size that can never exceed the population.
nreg <- max(1L, min(1000L, n_islands))
ntim <- 100L
if (n_islands < 10L)
  message("WARNING: very few islands (", n_islands, "); param estimates unstable")

ex <- tryCatch(
  ExoData(file = bam, mc.cores = mc_cores,
          nregions = nreg, ntimes = ntim, verbose = TRUE),
  error = function(e) {
    message("ExoData param pass failed (", conditionMessage(e),
            "); falling back to ex0"); ex0
  }
)

save_plot <- function(expr, path, w = 6, h = 5) {
  pdf(path, width = w, height = h)
  tryCatch(print(expr), error = function(e) {
    plot.new(); title(main = paste("plot unavailable:", conditionMessage(e)))
  })
  dev.off()
}

save_plot(ARCvURCplot(ex, names.input = sample_id), out_arc)
save_plot(FSRDistplot(ex, names.input = sample_id), out_fsr)
save_plot(regionCompplot(ex, names.input = sample_id), out_rc)
save_plot(paramDistBoxplot(ex, names.input = sample_id), out_pd)

## summary metrics (guarded)
b1 <- tryCatch(beta1(ex), error = function(e) NA_real_)
b2 <- tryCatch(beta2(ex), error = function(e) NA_real_)
med <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

fsr_median <- tryCatch(median(S4Vectors::mcols(ex)$FSR, na.rm = TRUE),
                       error = function(e) NA_real_)
nr <- tryCatch(sum(S4Vectors::mcols(ex)$depth), error = function(e) NA_real_)

stats <- data.frame(
  sample = sample_id, strain = strain, kind = kind,
  n_reads = nr, beta1_median = med(b1), beta2_median = med(b2),
  FSR_median = fsr_median, n_islands = n_islands,
  stringsAsFactors = FALSE
)
write.table(stats, file = out_stats, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote metrics: ", out_stats)
sink(type = "message"); sink(type = "output"); close(log)