## ============================================================================
## run_chipexoqual.R — ChIP-exo QC for one (pre-subsampled) BAM   [v5]
## ============================================================================
## Coefficient meaning — verified against BOTH the ChIPexoQual parameter
## documentation and the Welch et al. (NAR 2017) in-text definitions, which
## agree:
##   * beta1 = ENRICHMENT
##       "average number of aligned reads per unique position when sequencing
##        depth is large" -> high when reads pile onto few positions (a
##        strongly bound site), low when spread thinly (diffuse background).
##   * beta2 = LIBRARY COMPLEXITY
##       "overall change in depth as the width varies across regions"
##        (beta2 = -kappa/gamma) -> the axis separating complex libraries from
##        low-complexity / over-amplified ones. HIGH beta2 is the paper's
##        headline low-complexity warning.
##   Read together, at MATCHED DEPTH, as shape descriptors of the URC-vs-ARC
##   cloud — not as standalone quality scores.
##
## ExoData() forms ALL islands regardless of `nregions`; `nregions` only sets
## how many islands are sampled for the depth ~ 0 + uniquePos + width
## regression (repeated `ntimes`). The historical crash is the unguarded
## sample(.N, nregions) when nregions > n_islands.
##
## Input BAM is pre-subsampled to a fixed read count (rule subsample_bam_for_qc)
## because deep coverage on a ~2 Mb genome merges the chromosome into a single
## island. A fixed read count keeps beta1/beta2 comparable across samples.
##
## Outputs: ARC-URC, FSR, region-composition, paramDist boxplot, a DEPTH-LADDER
## panel (paper's signature diagnostic), and a one-row metrics TSV.
## ============================================================================

suppressMessages({
  library(ChIPexoQual)
  library(GenomicAlignments)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message"); sink(log, type = "output")

## ---- inputs / params / outputs --------------------------------------------
bam        <- snakemake@input[["bam"]]        # already subsampled to target depth
full_bam   <- snakemake@input[["full_bam"]]   # original (for the depth ladder)
sample_id  <- snakemake@params[["sample_id"]]
kind       <- snakemake@params[["kind"]]
strain     <- snakemake@params[["strain"]]
mc_cores   <- as.integer(snakemake@params[["mc_cores"]])
ladder     <- as.numeric(snakemake@params[["ladder_depths"]])  # e.g. c(1e5,2.5e5,5e5)

out_arc    <- snakemake@output[["arc_urc"]]
out_fsr    <- snakemake@output[["fsr"]]
out_rc     <- snakemake@output[["regioncomp"]]
out_pd     <- snakemake@output[["paramdist"]]
out_ladder <- snakemake@output[["ladder"]]
out_stats  <- snakemake@output[["stats"]]

message("ChIPexoQual on (subsampled): ", bam)

save_plot <- function(expr, path, w = 6, h = 5) {
  pdf(path, width = w, height = h)
  tryCatch(print(expr), error = function(e) {
    plot.new(); title(main = paste("plot unavailable:", conditionMessage(e)))
  })
  dev.off()
}

## ---- build the ExoData object ----------------------------------------------
## Cheap count pass (nregions=1 can never exceed the population) to size the
## regression safely, then the real build carrying the regression.
ex0 <- ExoData(file = bam, mc.cores = mc_cores,
               nregions = 1L, ntimes = 1L, verbose = FALSE)
n_islands <- length(ex0)
message("islands (subsampled bam): ", n_islands)

nreg <- min(1000L, max(2L, n_islands))
ex <- if (n_islands >= 2L) {
  ExoData(file = bam, mc.cores = mc_cores,
          nregions = nreg, ntimes = 100L, verbose = TRUE)
} else {
  ex0
}

## ---- island-level diagnostic panels ----------------------------------------
save_plot(ARCvURCplot(ex,     names.input = sample_id), out_arc)
save_plot(FSRDistplot(ex,     names.input = sample_id), out_fsr)
save_plot(regionCompplot(ex,  names.input = sample_id), out_rc)

## ---- beta1 / beta2 via the documented public accessors ---------------------
## beta1 = ENRICHMENT ; beta2 = LIBRARY COMPLEXITY
b1 <- tryCatch(beta1(ex), error = function(e) NA_real_)
b2 <- tryCatch(beta2(ex), error = function(e) NA_real_)

pd_tab <- tryCatch(paramDist(ex), error = function(e) NULL)
if (!is.null(pd_tab) && all(c("beta1", "beta2") %in% colnames(pd_tab))) {
  pdf(out_pd, width = 6, height = 5)
  boxplot(pd_tab$beta1, pd_tab$beta2,
          names = c("beta1\n(enrichment)", "beta2\n(complexity)"),
          main = sample_id, ylab = "coefficient estimate",
          col = c("#0E6E6B", "#A8471C"))
  dev.off()
} else {
  pdf(out_pd); plot.new(); title("paramDist unavailable"); dev.off()
}

## ---- DEPTH LADDER (paper's signature diagnostic) ---------------------------
## Built from the FULL bam so the ladder reflects true depth behaviour. Depths
## (config) stay in the regime where all samples form > nregions islands.
## Rising beta2 across increasing depth is the low-complexity warning.
ladder_ok <- FALSE
el <- tryCatch(
  ExoDataSubsampling(file = full_bam, sample.depth = sort(ladder),
                     nregions = 200L, ntimes = 100L,
                     mc.cores = mc_cores, verbose = FALSE),
  error = function(e) { message("ladder failed: ", conditionMessage(e)); NULL }
)
if (!is.null(el)) {
  save_plot(paramDistBoxplot(el), out_ladder)   # beta1/beta2 vs depth
  ladder_ok <- TRUE
} else {
  pdf(out_ladder); plot.new(); title("depth ladder unavailable"); dev.off()
}

## ---- summary metrics (single canonical depth for the scorecard) ------------
fsr_median <- tryCatch(median(S4Vectors::mcols(ex)$FSR, na.rm = TRUE),
                       error = function(e) NA_real_)
nr <- tryCatch(nreads(ex), error = function(e)
        tryCatch(sum(S4Vectors::mcols(ex)$depth), error = function(e2) NA_real_))

stats <- data.frame(
  sample = sample_id, strain = strain, kind = kind,
  n_reads_subsampled = nr,
  beta1_enrichment = b1, beta2_complexity = b2,
  FSR_median = fsr_median, n_islands = n_islands,
  ladder_generated = ladder_ok,
  stringsAsFactors = FALSE
)
write.table(stats, file = out_stats, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote metrics: ", out_stats)

sink(type = "message"); sink(type = "output"); close(log)