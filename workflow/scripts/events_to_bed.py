"""ChExMix events -> fixed-width BED windows, optionally filtered to a
high-confidence subset by replicate reproducibility.

peakset:
  "all"    -> every called binding event
  "hiconf" -> only events whose ChExMix replication code is a number >= min_reps
              (i.e. called significant in at least that many replicates).
              Codes 'X' (not called), 'S' (single replicate) and 'D'
              (significant difference between replicates) are excluded.

The join between the .events file and the .replicationcodes.table is done on
the event COORDINATE (not row order) — the two files are not guaranteed to be
in the same order.
"""
import os, re

events   = snakemake.input.events            # noqa: F821
repcodes = snakemake.input.repcodes          # noqa: F821
fai      = snakemake.input.fai               # noqa: F821
out      = snakemake.output.bed              # noqa: F821
half     = int(snakemake.params.half)        # noqa: F821
peakset  = str(snakemake.params.peakset)     # noqa: F821
min_reps = int(snakemake.params.min_reps)    # noqa: F821
logp     = snakemake.log[0]                  # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

# chromosome lengths
chrlen = {}
with open(fai) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            chrlen[f[0]] = int(f[1])

# replication code per event coordinate
code = {}
with open(repcodes) as fh:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip() != "X":
            code[f[0].strip()] = f[1].strip()

def passes(point):
    if peakset == "all":
        return True
    c = code.get(point)
    if c is None:
        return False
    try:
        return int(c) >= min_reps       # numeric codes = n replicates
    except ValueError:
        return False                     # 'S' or 'D' excluded from hiconf

n_in = n_out = n_nojoin = 0
with open(events) as fh, open(out, "w") as ob, open(logp, "w") as lg:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        n_in += 1
        point = ln.split("\t")[0].strip()
        if point not in code:
            n_nojoin += 1
        if not passes(point):
            continue
        m = re.match(r"^(?:chr)?(.+):(\d+)$", point)
        if not m:
            lg.write(f"skip (unparsed): {point}\n"); continue
        chrom, pos = m.group(1), int(m.group(2))
        if chrom not in chrlen and ("chr" + chrom) in chrlen:
            chrom = "chr" + chrom
        if chrom not in chrlen:
            lg.write(f"skip (chrom not in fai): {chrom}\n"); continue
        start = max(0, (pos - 1) - half)
        end   = min(chrlen[chrom], (pos - 1) + half)
        ob.write(f"{chrom}\t{start}\t{end}\t{chrom}_{pos}\t.\t+\n")
        n_out += 1
    lg.write(f"peakset={peakset} min_reps={min_reps} events_in={n_in} "
             f"windows_out={n_out} unjoined_coords={n_nojoin} half={half}\n")