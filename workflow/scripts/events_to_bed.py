"""Convert ChExMix .events Point coordinates to fixed-width BED windows.

ChExMix event Point format (from the .events file, non-# lines):
    chrNC_007432.1:874875   ...
i.e. 'chr' + <chromName> + ':' + <1-based position>, tab-separated columns.
We emit BED (0-based, half-open) windows of ±half around each Point, clipped
to chromosome ends read from the .fai.
"""
import os, re

events = snakemake.input.events          # noqa: F821
fai    = snakemake.input.fai             # noqa: F821
out    = snakemake.output.bed            # noqa: F821
half   = int(snakemake.params.half)      # noqa: F821
logp   = snakemake.log[0]                # noqa: F821

# chromosome lengths from .fai (col1 name, col2 length)
chrlen = {}
with open(fai) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            chrlen[f[0]] = int(f[1])

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

n_in = n_out = 0
with open(events) as fh, open(out, "w") as ob, open(logp, "w") as lg:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        n_in += 1
        point = ln.split("\t")[0].strip()   # e.g. chrNC_007432.1:874875
        m = re.match(r"^(?:chr)?(.+):(\d+)$", point)
        if not m:
            lg.write(f"skip (unparsed): {point}\n"); continue
        chrom, pos = m.group(1), int(m.group(2))
        # strip a leading 'chr' if the fai names lack it; try both
        if chrom not in chrlen and ("chr" + chrom) in chrlen:
            chrom = "chr" + chrom
        if chrom not in chrlen:
            # ChExMix prefixes 'chr'; fai uses bare accession -> already stripped above.
            lg.write(f"skip (chrom not in fai): {chrom}\n"); continue
        start = max(0, (pos - 1) - half)        # pos is 1-based -> 0-based center
        end   = min(chrlen[chrom], (pos - 1) + half)
        name  = f"{chrom}_{pos}"
        ob.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t+\n")
        n_out += 1
    lg.write(f"events_in={n_in} windows_out={n_out} half={half}\n")