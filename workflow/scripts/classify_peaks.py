"""Classify each ChExMix peak as canonical (GAS Per box present in its window)
or non-canonical, and attach the ChExMix subtype.

Inputs:
  windows : BED of peak windows (name = chrom_pos, matches event Point)
  events  : ChExMix .events (Point, ... SubtypeName column)
  fimo    : FIMO tsv of GAS Per-box hits genome-wide
Output:
  per-peak TSV: chrom, point, subtype, perbox_in_window(bool), n_perbox_hits,
                best_fimo_qvalue
"""
import os, re, csv

windows = snakemake.input.windows        # noqa: F821
events  = snakemake.input.events         # noqa: F821
fimo    = snakemake.input.fimo           # noqa: F821
out     = snakemake.output.tsv           # noqa: F821
logp    = snakemake.log[0]               # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

# 1. windows: name -> (chrom, start, end)
win = {}
with open(windows) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4: continue
        win[f[3]] = (f[0], int(f[1]), int(f[2]))

# 2. events: Point -> subtype (SubtypeName column, if present)
#    header line (starts with #Point) tells us the column index.
subtype = {}
sub_idx = None
with open(events) as fh:
    for ln in fh:
        if ln.startswith("#Point"):
            cols = ln.lstrip("#").rstrip("\n").split("\t")
            for i, c in enumerate(cols):
                if c.strip().lower() == "subtypename":
                    sub_idx = i
            continue
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        point = f[0].strip()
        m = re.match(r"^(?:chr)?(.+):(\d+)$", point)
        if not m: continue
        key = f"{('chr'+m.group(1)) if not m.group(1).startswith('chr') else m.group(1)}_{m.group(2)}"
        # match the naming used in events_to_bed (chrom kept as-is there)
        key_bare = f"{m.group(1)}_{m.group(2)}"
        st = f[sub_idx].strip() if (sub_idx is not None and len(f) > sub_idx) else "NA"
        subtype[key] = st
        subtype[key_bare] = st

# 3. FIMO hits: list of (seqname, start, stop, qvalue)
fimo_hits = []
with open(fimo) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for row in rdr:
        if not row or row.get("motif_id", "").startswith("#"): continue
        seq = row.get("sequence_name") or row.get("sequence name")
        try:
            start = int(row["start"]); stop = int(row["stop"])
            qv = float(row.get("q-value", row.get("q_value", "nan")))
        except (KeyError, ValueError, TypeError):
            continue
        fimo_hits.append((seq, start, stop, qv))

def hits_in(chrom, s, e):
    got = []
    for (seq, hs, he, qv) in fimo_hits:
        # FIMO is 1-based inclusive; window is 0-based half-open. Overlap test.
        if seq == chrom and not (he < s + 1 or hs > e):
            got.append(qv)
    return got

with open(out, "w") as ob, open(logp, "w") as lg:
    ob.write("peak\tchrom\twindow_start\twindow_end\tsubtype\t"
             "perbox_in_window\tn_perbox_hits\tbest_qvalue\n")
    n_canon = 0
    for name, (chrom, s, e) in sorted(win.items()):
        qs = hits_in(chrom, s, e)
        canon = len(qs) > 0
        n_canon += int(canon)
        st = subtype.get(name, "NA")
        best = min(qs) if qs else ""
        ob.write(f"{name}\t{chrom}\t{s}\t{e}\t{st}\t{canon}\t{len(qs)}\t{best}\n")
    lg.write(f"peaks={len(win)} canonical={n_canon} "
             f"fimo_hits_total={len(fimo_hits)}\n")