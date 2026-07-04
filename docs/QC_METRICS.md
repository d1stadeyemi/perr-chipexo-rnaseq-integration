# QC metrics — what each number means and why

This stage certifies the sequencing-core delivery **before** any peak calling.
Every threshold below is advisory and sourced; ChIP-exo deliberately violates
several ChIP-seq QC conventions, so we report values and gentle flags rather
than hard pass/fail on complexity.

## 1. Genome-build verification (`header_check`)
The BAM `@SQ` header records the exact reference the core aligned to. We assert
it against the accession declared in `config.yaml`:

| strain | serotype | accession |
|---|---|---|
| A909 | Ia | `NC_007432.1` (2,127,839 bp, complete) |
| CJB111 | V | `NZ_CP063198.2` (2021 closed assembly) |

A mismatch **aborts** the run — using the wrong build silently is the classic
way to invalidate a ChIP-exo analysis. The `@PG` line also records the aligner
and whether Picard/MarkDuplicates was applied (the delivery is pre-deduplicated).

## 2. Mapping & pairing (`samtools flagstat`)
Reported per sample. The modENCODE point-source objective is ≥2 M uniquely
mapped reads per replicate (Landt et al., *Genome Res.* 2012; ENCODE/modENCODE
guidelines) — all eight samples exceed this by an order of magnitude. Lower
mapping rates in CJB111 vs A909 are expected (greater divergence from its
reference and a half-scale IP; see the core report) and are informational, not
disqualifying.

## 3. FastQC
Raw-read quality and adapter content. Given the short insert sizes noted in the
core report (~67–85 bp average), R1/R2 overlap heavily and some adapter
read-through is expected — this is diagnostic context if you ever re-process
from FASTQ, not a failure of the provided BAMs.

## 4. ChIPexoQual (the ChIP-exo-specific QC) — Welch et al., *NAR* 2017
Computed per BAM over "read islands" (enriched regions):

- **ARC vs URC plot** — Average Read Coefficient vs Unique Read Coefficient.
  A good ChIP-exo library shows the characteristic two-arm structure; the
  arms encode the balance of enrichment and complexity.
- **FSR (Forward Strand Ratio) distribution** — strand balance. Values near
  0.5 (config window 0.4–0.6) indicate the forward/reverse symmetry that
  strand-aware ChIP-exo callers (GEM, MACE, CexoR) depend on. Strong imbalance
  is flagged `STRAND_IMBALANCE`.
- **beta1** — regression coefficient tracking **enrichment** (higher = better).
- **beta2** — coefficient tracking **library complexity**; larger beta2 flags a
  low-complexity library (Welch et al. show beta2 rising for their lower-quality
  replicates even at high depth). Reported, not hard-failed.

## 5. deepTools grouped QC (per strain: 3 IP reps + 1 control)
- **plotFingerprint** — cumulative-enrichment curve. A specific IP is strongly
  skewed away from the diagonal; the no-tag control sits near the diagonal.
  The IP-vs-control separation is the single most intuitive "did the pulldown
  work" plot and mirrors the tagged-vs-untagged contrast in the core report.
- **multiBamSummary + plotCorrelation** (Spearman) — IP replicates should
  correlate tightly and cluster apart from the control.
- **bamPEFragmentSize** — the paired-end insert-size distribution, which should
  reproduce the short ChIP-exo fragments the core observed.

## On deduplication (important caveat)
The core provides **deduplicated** BAMs. Carroll et al. (2014) showed naive
PCR-duplicate removal in ChIP-exo can strip real signal, because reads pile up
at identical 5′ positions by biology, not only by PCR. Paired-end full-fragment
deduplication (matching both ends) mitigates this. Because the FASTQs are
retained, a later optional module re-derives peaks with and without dedup as a
**sensitivity analysis** — turning this caveat into a reported robustness check
rather than an unstated assumption.

## Design note: 3 IP replicates, 1 control
IP is triplicated to establish reproducibility of binding (where the variance
lives); the no-tag control models background, which is far less variable, so a
single control per strain/condition is standard and explicitly permitted by
ENCODE. The tradeoff — no estimate of control-to-control variability — is a
one-line limitation.
