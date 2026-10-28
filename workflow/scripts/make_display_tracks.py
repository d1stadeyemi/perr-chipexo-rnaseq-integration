#!/usr/bin/env python3
"""Display-only IGV tracks for ChIP-exo figures.

Writes two kinds of track from existing outputs. Neither is used in any
analysis; both are for figures, and the legend must say what they are.

1. EVENT S/N TRACK (recommended; matches the "S/N" rows in Seo et al. 2015,
   Nat Commun, GadE/W/X).
   One bar per ChExMix event, height = IP/control fold enrichment, width
   2*halfwin+1 bp. Zero everywhere else. This is a picture of the peak calls,
   not of the reads.

2. MASKED 5'-END TRACKS
   The strand-specific 5'-end bigwigs with every base farther than `halfwin`
   from a called event set to zero. Read heights inside the windows are
   unchanged. The hidden reads still exist, so the legend must say
   "signal shown only within +/-N bp of called events".

Usage (run inside the pileup/deeptools env, which provides pyBigWig):

  python make_display_tracks.py \
      --events results/03_peaks/chexmix/A909/A909.events \
      --bw results/10_igv/bigwig/<depth>/A909/A909.ip.fwd.cpm.bw \\
           results/10_igv/bigwig/<depth>/A909/A909.ip.rev.neg.cpm.bw \
      --halfwin 50 --out igv_display/A909

Outputs (bedGraph; IGV loads these directly, no conversion needed):
  <out>.events_SN.bedgraph
  <out>.<bigwig basename>.masked.bedgraph   (one per --bw)
"""
import argparse
import os
import re
import sys

import pyBigWig


def read_events(path):
    """Return [(chrom, pos_1based, fold)] from a ChExMix .events file.

    ChExMix files open with '###' comment lines; the column header is the line
    starting '#Point'. The fold column is the first header field containing
    'log2fold'; if none is found, column 4 (index 3) is used, which is where
    make_igv_tracks.py already reads log2FC. A leading 'chr' on the chromosome
    is stripped so names match the bigwigs (NC_007432.1)."""
    out, fold_i = [], 3
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#Point"):
                cols = ln.lstrip("#").rstrip("\n").split("\t")
                hits = [i for i, c in enumerate(cols) if "log2fold" in c.lower()]
                if hits:
                    fold_i = hits[0]
                continue
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0].strip())
            if not m:
                continue
            try:
                fold = 2 ** float(f[fold_i])
            except (IndexError, ValueError):
                continue
            out.append((m.group(1), int(m.group(2)), fold))
    if not out:
        sys.exit(f"No events parsed from {path}")
    return out


def merged_windows(events, halfwin, sizes):
    """Event windows as 0-based half-open intervals, merged per chromosome."""
    by = {}
    for c, p, _ in events:
        s = max(0, p - 1 - halfwin)
        e = min(sizes[c], p + halfwin)
        by.setdefault(c, []).append((s, e))
    merged = {}
    for c, iv in by.items():
        iv.sort()
        cur = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= cur[-1][1]:
                cur[-1][1] = max(cur[-1][1], e)
            else:
                cur.append([s, e])
        merged[c] = cur
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--bw", nargs="+", required=True)
    ap.add_argument("--halfwin", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    events = read_events(a.events)
    sizes = pyBigWig.open(a.bw[0]).chroms()
    missing = {c for c, _, _ in events} - set(sizes)
    if missing:
        sys.exit(f"Event chromosomes {sorted(missing)} not in bigwig "
                 f"chromosomes {sorted(sizes)}")
    wins = merged_windows(events, a.halfwin, sizes)

    # 1. event S/N track
    sn = f"{a.out}.events_SN.bedgraph"
    with open(sn, "w") as o:
        o.write(f'track type=bedGraph name="events S/N (fold, +/-{a.halfwin}bp)"\n')
        for c, p, fold in sorted(events):
            s, e = max(0, p - 1 - a.halfwin), min(sizes[c], p + a.halfwin)
            o.write(f"{c}\t{s}\t{e}\t{fold:.3f}\n")
    print(f"{len(events)} events -> {sn}")

    # 2. masked 5'-end tracks
    for path in a.bw:
        bw = pyBigWig.open(path)
        name = os.path.basename(path).replace(".bw", "").replace(".bigwig", "")
        dst = f"{a.out}.{name}.masked.bedgraph"
        kept = total = 0.0
        for c, L in bw.chroms().items():
            total += sum(abs(v) * (e - s) for s, e, v in (bw.intervals(c) or []))
        with open(dst, "w") as o:
            o.write(f'track type=bedGraph name="{name} (+/-{a.halfwin}bp of events)"\n')
            for c in sorted(wins):
                for s, e in wins[c]:
                    for is_, ie, v in bw.intervals(c, s, e) or []:
                        is_, ie = max(is_, s), min(ie, e)
                        if v != 0 and ie > is_:
                            o.write(f"{c}\t{is_}\t{ie}\t{v:.4f}\n")
                            kept += abs(v) * (ie - is_)
        bw.close()
        frac = kept / total if total else 0
        print(f"{name}: {frac:.1%} of 5'-end signal lies inside event windows -> {dst}")


if __name__ == "__main__":
    main()
