"""Classify each ChExMix peak by Per-box content, and attach the ChExMix subtype.

TWO MOTIFS ARE SCORED INDEPENDENTLY:

  GAS box (TTANAATNATTNTAA)  the Streptococcus pyogenes consensus. A REFERENCE
                            YARDSTICK ONLY. This project showed the GBS motif
                            is NOT the GAS box (TOMTOM p = 0.34, 6 bp overlap),
                            so GAS-box hits are reported for comparison and are
                            not evidence of PerR binding in GBS.

  GBS motif (TTMTAAATAAG)   discovered DE NOVO under the peaks in this project
                            and recovered independently by MEME and STREME in
                            both serotypes. This is the biologically relevant
                            motif and the one that should drive classification.

WHY A MERGED IUPAC CONSENSUS RATHER THAN EACH STRAIN'S OWN PWM: MEME reported
A909 on the minus strand (CTTATTTAKAA, the reverse complement) and CJB111 on
the plus strand (TTCTAAATAAS); STREME recovered a long degenerate form in the
A909 hiconf set (AKWAWAAATAAYHNTW) and a 12-mer in CJB111 (ATTCTAAATAAG). All
four agree on the core TAAATAA and differ only at the terminal position.
Scanning each genome with its own matrix would compare counts obtained with two
DIFFERENT matrices; one consensus scanned identically in both genomes is the
only clean cross-strain comparison. The strain-specific PWMs remain available
in memechip_*/combined.meme as a consistency check.

Inputs:
  windows  : BED of peak windows (name = chrom_pos, matches the event Point)
  events   : ChExMix .events (Point, ... SubtypeName column)
  fimo     : FIMO tsv, GAS box, genome-wide
  fimo_gbs : FIMO tsv, GBS motif, genome-wide
Output, one row per peak:
  peak, chrom, window_start, window_end, subtype,
  perbox_in_window, n_perbox_hits, best_qvalue,                    <- GAS
  gbs_perbox_in_window, n_gbs_perbox_hits, best_gbs_qvalue,        <- GBS
  gbs_best_offset                                                  <- GBS

`gbs_best_offset` is the signed distance in bp from the ChExMix event Point to
the CENTRE of the best-scoring GBS motif hit in the window (negative = motif
lies 5' of the crosslink on the plus strand of the reference). For a genuine
ChIP-exo footprint this should cluster tightly around a small fixed offset
rather than scattering across the window, which is the positional-enrichment
statistic CentriMo would have produced had it not crashed in this MEME build.

NOTE: the GAS columns keep their v1 names and positions so that any downstream
positional indexing continues to work; the four GBS columns are APPENDED.
"""
import os
import re
import csv

windows = snakemake.input.windows        # noqa: F821
events = snakemake.input.events          # noqa: F821
fimo = snakemake.input.fimo              # noqa: F821
fimo_gbs = snakemake.input.fimo_gbs      # noqa: F821
out = snakemake.output.tsv               # noqa: F821
logp = snakemake.log[0]                  # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. windows: name -> (chrom, start, end)
# ---------------------------------------------------------------------------
win = {}
with open(windows) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        win[f[3]] = (f[0], int(f[1]), int(f[2]))

# ---------------------------------------------------------------------------
# 2. events: Point -> subtype (SubtypeName column, located from the #Point line)
# ---------------------------------------------------------------------------
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
        if not m:
            continue
        base = m.group(1)
        key_chr = f"{('chr' + base) if not base.startswith('chr') else base}_{m.group(2)}"
        key_bare = f"{base}_{m.group(2)}"
        st = f[sub_idx].strip() if (sub_idx is not None and len(f) > sub_idx) else "NA"
        subtype[key_chr] = st
        subtype[key_bare] = st

warnings = []
if sub_idx is None:
    warnings.append(
        "WARNING: no '#Point' header line found in the events file; "
        "SubtypeName column could not be located and every peak will "
        "be reported as subtype=NA")


# ---------------------------------------------------------------------------
# 3. FIMO hit tables
# ---------------------------------------------------------------------------
def read_fimo(path):
    """-> list of (seqname, start, stop, qvalue). FIMO is 1-based inclusive."""
    hits = []
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for rowd in rdr:
            if not rowd or str(rowd.get("motif_id", "")).startswith("#"):
                continue
            seq = rowd.get("sequence_name") or rowd.get("sequence name")
            try:
                start = int(rowd["start"])
                stop = int(rowd["stop"])
                qv = float(rowd.get("q-value", rowd.get("q_value", "nan")))
            except (KeyError, ValueError, TypeError):
                continue
            hits.append((seq, start, stop, qv))
    return hits


gas_hits = read_fimo(fimo)
gbs_hits = read_fimo(fimo_gbs)


def hits_in(hit_list, chrom, s, e):
    """Overlap between a 1-based inclusive FIMO hit and a 0-based half-open
    BED window. Returns [(qvalue, hit_centre_1based), ...]."""
    got = []
    for (seq, hs, he, qv) in hit_list:
        if seq == chrom and not (he < s + 1 or hs > e):
            got.append((qv, (hs + he) / 2.0))
    return got


# ---------------------------------------------------------------------------
# 4. write
# ---------------------------------------------------------------------------
HDR = ["peak", "chrom", "window_start", "window_end", "subtype",
       "perbox_in_window", "n_perbox_hits", "best_qvalue",
       "gbs_perbox_in_window", "n_gbs_perbox_hits", "best_gbs_qvalue",
       "gbs_best_offset"]

with open(out, "w") as ob, open(logp, "w") as lg:
    for w in warnings:
        lg.write(w + "\n")
    ob.write("\t".join(HDR) + "\n")
    n_gas = n_gbs = n_both = 0
    offsets = []

    for name, (chrom, s, e) in sorted(win.items()):
        # the window name is chrom_point; recover the ChExMix Point
        try:
            point = int(name.rsplit("_", 1)[1])
        except (IndexError, ValueError):
            point = (s + e) // 2          # fall back to window centre

        g_hits = hits_in(gas_hits, chrom, s, e)
        b_hits = hits_in(gbs_hits, chrom, s, e)
        g_can, b_can = len(g_hits) > 0, len(b_hits) > 0
        n_gas += int(g_can)
        n_gbs += int(b_can)
        n_both += int(g_can and b_can)

        g_best = min(q for q, _ in g_hits) if g_hits else ""
        if b_hits:
            bq, bc = min(b_hits, key=lambda x: x[0])
            b_best = bq
            b_off = int(round(bc - point))
            offsets.append(b_off)
        else:
            b_best, b_off = "", ""

        ob.write("\t".join(str(x) for x in [
            name, chrom, s, e, subtype.get(name, "NA"),
            g_can, len(g_hits), g_best,
            b_can, len(b_hits), b_best, b_off]) + "\n")

    lg.write(f"peaks={len(win)}\n")
    lg.write(f"GAS-box hits total={len(gas_hits)}  peaks with a GAS box={n_gas}\n")
    lg.write(f"GBS-motif hits total={len(gbs_hits)}  peaks with a GBS motif={n_gbs}\n")
    lg.write(f"peaks with BOTH={n_both}\n")
    if offsets:
        offsets.sort()
        mid = offsets[len(offsets) // 2]
        within50 = sum(1 for o in offsets if abs(o) <= 50)
        lg.write(f"GBS motif offset from crosslink: median={mid} bp, "
                 f"within +/-50 bp={within50}/{len(offsets)}\n")