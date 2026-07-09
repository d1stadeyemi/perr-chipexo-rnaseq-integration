"""Collate QC outputs into a single decision scorecard.

Combines, per sample:
  * genome-build verification  (check_bam_header)
  * mapping / pairing          (samtools flagstat)
  * ChIP-exo QC                (ChIPexoQual: beta1, beta2, FSR)

Flags are ADVISORY. ChIP-exo deliberately breaks several ChIP-seq QC rules of
thumb (it is expected to be highly duplicated at footprint borders), so we
report values and gentle flags rather than hard pass/fail on complexity.
See docs/QC_METRICS.md for the rationale and citations.
"""

import os
import re

samples = snakemake.params.samples                    # noqa: F821
meta = snakemake.params.meta                           # noqa: F821
min_um = int(snakemake.params.min_uniquely_mapped)     # noqa: F821
fsr_low = float(snakemake.params.fsr_low)              # noqa: F821
fsr_high = float(snakemake.params.fsr_high)            # noqa: F821

out_tsv = snakemake.output.tsv                          # noqa: F821
out_md = snakemake.output.md                            # noqa: F821
logpath = snakemake.log[0]                              # noqa: F821

qc_dir = os.path.dirname(os.path.dirname(snakemake.input.chipexoqual[0]))  # noqa: F821
header_dir = os.path.dirname(snakemake.input.header[0])                     # noqa: F821
samtools_dir = os.path.dirname(snakemake.input.flagstat[0])                 # noqa: F821

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
os.makedirs(os.path.dirname(logpath), exist_ok=True)


def read_flagstat(path):
    """Return dict with total, mapped, mapped_pct, properly_paired_pct, dup."""
    d = {"mapped": None, "mapped_pct": None, "properly_paired_pct": None,
         "primary_dup": None, "total": None}
    with open(path) as fh:
        for ln in fh:
            if "in total" in ln:
                d["total"] = int(ln.split()[0])
            elif re.search(r"\bmapped\b", ln) and "primary" not in ln and "mate" not in ln:
                d["mapped"] = int(ln.split()[0])
                m = re.search(r"\(([\d.]+)%", ln)
                if m:
                    d["mapped_pct"] = float(m.group(1))
            elif "properly paired" in ln:
                m = re.search(r"\(([\d.]+)%", ln)
                if m:
                    d["properly_paired_pct"] = float(m.group(1))
            elif "duplicates" in ln and "primary duplicates" not in ln:
                d["primary_dup"] = int(ln.split()[0])
    return d


def read_header(path):
    with open(path) as fh:
        fh.readline()  # header row
        row = fh.readline().rstrip("\n").split("\t")
    # sample, expected_acc, observed_SN, observed_LN, aligner, status, notes
    return {
        "expected_accession": row[1],
        "observed_SN": row[2],
        "aligner": row[4],
        "build_status": row[5],
    }


def read_chipexoqual(path):
    with open(path) as fh:
        keys = fh.readline().rstrip("\n").split("\t")
        vals = fh.readline().rstrip("\n").split("\t")
    return dict(zip(keys, vals))


rows = []
for s in samples:
    m = meta[s]
    fs = read_flagstat(os.path.join(samtools_dir, f"{s}.flagstat.txt"))
    hd = read_header(os.path.join(header_dir, f"{s}.header_check.tsv"))
    cq = read_chipexoqual(os.path.join(qc_dir, s, f"{s}.metrics.tsv"))

    def fnum(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None

    fsr = fnum(cq.get("FSR_median"))
    beta1 = fnum(cq.get("beta1_median"))
    beta2 = fnum(cq.get("beta2_median"))
    nreads = fnum(cq.get("n_reads_subsampled"))

    flags = []
    if hd["build_status"] != "PASS":
        flags.append("GENOME_BUILD_MISMATCH")
    if fs["mapped"] is not None and fs["mapped"] < min_um:
        flags.append("LOW_MAPPED")
    if fsr is not None and not (fsr_low <= fsr <= fsr_high):
        flags.append("STRAND_IMBALANCE")

    rows.append({
        "sample": s,
        "strain": m["strain"],
        "kind": m["kind"],
        "replicate": m["replicate"],
        "accession": m["accession"],
        "build_status": hd["build_status"],
        "aligner": hd["aligner"],
        "total_reads": fs["total"],
        "mapped": fs["mapped"],
        "mapped_pct": fs["mapped_pct"],
        "properly_paired_pct": fs["properly_paired_pct"],
        "chipexo_nreads": nreads,
        "beta1_enrichment": beta1,
        "beta2_lowcomplexity": beta2,
        "FSR_median": fsr,
        "flags": ";".join(flags) if flags else "ok",
    })

cols = ["sample", "strain", "kind", "replicate", "accession", "build_status",
        "aligner", "total_reads", "mapped", "mapped_pct",
        "properly_paired_pct", "chipexo_nreads", "beta1_enrichment",
        "beta2_lowcomplexity", "FSR_median", "flags"]

with open(out_tsv, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

# --- readable markdown --------------------------------------------------------
with open(out_md, "w") as fh:
    fh.write("# ChIP-exo QC scorecard\n\n")
    fh.write("Flags are advisory. See `docs/QC_METRICS.md` for thresholds and "
             "why ChIP-exo tolerates high duplication.\n\n")
    fh.write("| sample | strain | kind | build | mapped % | prop.paired % | "
             "beta1 (enrich) | beta2 (low-cplx) | FSR | flags |\n")
    fh.write("|---|---|---|---|---|---|---|---|---|---|\n")
    for r in rows:
        fh.write(
            f"| {r['sample']} | {r['strain']} | {r['kind']} | "
            f"{r['build_status']} | {r['mapped_pct']} | "
            f"{r['properly_paired_pct']} | {r['beta1_enrichment']} | "
            f"{r['beta2_lowcomplexity']} | {r['FSR_median']} | {r['flags']} |\n"
        )

with open(logpath, "w") as log:
    log.write(f"scorecard written for {len(rows)} samples\n")
