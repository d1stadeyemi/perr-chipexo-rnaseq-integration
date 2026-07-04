"""Verify the reference build recorded in the BAM header.

Reads @SQ (sequence name + length) and @PG (aligner / dedup provenance) from
the BAM header via `samtools view -H` and checks the sequence name against the
accession declared in config for this sample's strain. A mismatch aborts the
workflow — using the wrong genome build silently is the classic way to ruin a
ChIP-exo analysis.
"""

import os
import subprocess
import sys

bam = snakemake.input.bam                    # noqa: F821
out = snakemake.output.report                # noqa: F821
accession = snakemake.params.accession       # noqa: F821
chrom_len = snakemake.params.chrom_len       # noqa: F821
sample_id = snakemake.params.sample_id       # noqa: F821
logpath = snakemake.log[0]                    # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logpath), exist_ok=True)

header = subprocess.run(
    ["samtools", "view", "-H", bam],
    capture_output=True, text=True, check=True,
).stdout

sq_lines = [ln for ln in header.splitlines() if ln.startswith("@SQ")]
pg_lines = [ln for ln in header.splitlines() if ln.startswith("@PG")]

sq = {}
for ln in sq_lines:
    parts = dict(
        p.split(":", 1) for p in ln.split("\t")[1:] if ":" in p
    )
    if "SN" in parts:
        sq[parts["SN"]] = parts.get("LN")

# aligner provenance (bwa? bowtie2? Picard MarkDuplicates?)
programs = []
for ln in pg_lines:
    parts = dict(p.split(":", 1) for p in ln.split("\t")[1:] if ":" in p)
    programs.append(parts.get("PN", parts.get("ID", "?")) + " " + parts.get("VN", ""))

# --- checks -------------------------------------------------------------------
# Accept an exact SN match, or a match ignoring version suffix (NC_007432 vs .1)
def base_acc(x):
    return x.split(".")[0]

sn_names = list(sq.keys())
match = any(
    sn == accession or base_acc(sn) == base_acc(accession) for sn in sn_names
)

problems = []
if not match:
    problems.append(
        f"@SQ SN {sn_names} does not match expected accession '{accession}'"
    )
if chrom_len is not None and match:
    ln = None
    for sn in sn_names:
        if sn == accession or base_acc(sn) == base_acc(accession):
            ln = sq[sn]
    if ln is not None and str(ln) != str(chrom_len):
        problems.append(
            f"chromosome length {ln} != expected {chrom_len} (soft check)"
        )

status = "PASS" if not any("does not match" in p for p in problems) else "FAIL"

with open(out, "w") as fh:
    fh.write("sample\texpected_accession\tobserved_SN\tobserved_LN\t"
             "aligner_provenance\tstatus\tnotes\n")
    fh.write(
        f"{sample_id}\t{accession}\t{','.join(sn_names)}\t"
        f"{','.join(str(sq[s]) for s in sn_names)}\t"
        f"{'; '.join(p.strip() for p in programs)}\t{status}\t"
        f"{'; '.join(problems) if problems else 'ok'}\n"
    )

with open(logpath, "w") as log:
    log.write(header)

if status == "FAIL":
    sys.exit(
        f"[check_bam_header] {sample_id}: genome build mismatch — {problems}"
    )
