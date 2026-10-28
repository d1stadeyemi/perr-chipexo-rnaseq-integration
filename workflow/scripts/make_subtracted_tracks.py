#!/usr/bin/env python3
"""Control-subtracted ChIP-exo 5'-end tracks for IGV display.

For every base, on each strand:

    shown = max( IP  -  scale x background , 0 )

where `background` is the no-tag control at that base, averaged over a sliding
window of --smooth bp (default 500). Reverse-strand output is negated for the
mirrored display.

WHY THE CONTROL IS SMOOTHED
  At 1 bp, a control with "6 reads" at a base is itself a random draw. The
  background expected in the IP at that base is the control's local AVERAGE,
  not the exact count at that base: most bases where the IP has a read have
  no control read, so a literal base-by-base subtraction removes almost
  nothing. Smoothing gives the expected background, which is what MACS2 does
  with its local lambda. Use --smooth 1 to see the literal version.

WHY THE SCALE IS ESTIMATED
  The IP strand tracks count R1 5' ends on one strand; the control track counts
  all mates on both strands, so the two CPM tracks are not on the same per-base
  scale. `scale` = total IP / total control over bases farther than
  --exclude bp from any called event, so it is fitted on background only and
  calculated separately for each strand. It is printed so it can go in the legend.

WHAT THIS DOES NOT DO
  It removes the AVERAGE background, not random noise around that average.
  Bases where the IP happened to draw more reads than average stay above zero.
  The track is cleaner, not empty. Display only: not used by any analysis.

Usage (needs numpy + pyBigWig):
  python make_subtracted_tracks.py \
     --fwd  results/10_igv/bigwig/full_depth/A909/A909.ip.fwd.cpm.bw \
     --rev  results/10_igv/bigwig/full_depth/A909/A909.ip.rev.cpm.bw \
     --ctrl results/10_igv/bigwig/full_depth/A909/A909.ctrl.cpm.bw \
     --events results/03_peaks/chexmix/A909/A909.events \
     --out  results/igv_display_new/A909.sub
"""
import argparse
import os
import re

import numpy as np
import pyBigWig


def load(path):
    bw = pyBigWig.open(path)
    out = {c: np.nan_to_num(np.abs(np.array(bw.values(c, 0, L), dtype=float)))
           for c, L in bw.chroms().items()}
    bw.close()
    return out


def smooth(x, w):
    if w <= 1:
        return x
    k = np.ones(w) / w
    return np.convolve(x, k, mode="same")


def event_mask(path, sizes, half):
    """True at bases within `half` bp of any event (ChExMix #Point format)."""
    m = {c: np.zeros(L, bool) for c, L in sizes.items()}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            g = re.match(r"^(?:chr)?(.+):(\d+)", ln.split("\t")[0].strip())
            if g and g.group(1) in m:
                p = int(g.group(2)) - 1
                m[g.group(1)][max(0, p - half):p + half + 1] = True
    return m


def write_bedgraph(arrs, path, name, sign):
    with open(path, "w") as o:
        o.write(f'track type=bedGraph name="{name}"\n')
        for c, v in arrs.items():
            v = np.round(v, 4)
            nz = np.flatnonzero(v)
            if nz.size == 0:
                continue
            # run-length encode consecutive identical non-zero values
            brk = np.flatnonzero((np.diff(nz) != 1) | (np.diff(v[nz]) != 0)) + 1
            for seg in np.split(nz, brk):
                o.write(f"{c}\t{seg[0]}\t{seg[-1] + 1}\t{sign * v[seg[0]]:.4f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fwd", required=True)
    ap.add_argument("--rev", required=True, help="the NON-negated reverse track")
    ap.add_argument("--ctrl", required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--smooth", type=int, default=500)
    ap.add_argument("--exclude", type=int, default=500)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    ctrl = load(a.ctrl)
    sizes = {c: len(v) for c, v in ctrl.items()}
    near = event_mask(a.events, sizes, a.exclude)
    bg = {c: smooth(v, a.smooth) for c, v in ctrl.items()}

    for label, path, sign in (("fwd", a.fwd, 1), ("rev", a.rev, -1)):
        ip = load(path)
        far_ip = sum(ip[c][~near[c]].sum() for c in ip)
        far_bg = sum(bg[c][~near[c]].sum() for c in ip)
        scale = far_ip / far_bg if far_bg else 1.0
        sub = {c: np.clip(ip[c] - scale * bg[c], 0, None) for c in ip}
        kept = sum(s.sum() for s in sub.values()) / sum(v.sum() for v in ip.values())
        dst = f"{a.out}.{label}.bedgraph"
        write_bedgraph(sub, dst,
                       f"IP {label} 5' minus control (smooth {a.smooth}bp)", sign)
        print(f"{label}: scale={scale:.3f}  signal remaining={kept:.1%}  -> {dst}")


if __name__ == "__main__":
    main()
