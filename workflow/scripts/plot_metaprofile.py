"""plot_metaprofile.py -- the canonical ChIP-exo footprint figure.

WHAT IT SHOWS
  Average strand-separated 5' end signal, aligned on every ChExMix crosslink.
  Lambda-exonuclease digests until it collides with the crosslinked protein, so
  forward-strand 5' ends accumulate on the LEFT border of the footprint and
  reverse-strand ends on the RIGHT. The signature is therefore a forward peak
  just left of zero and a reverse peak just right of it, with the separation
  giving the footprint width.

  This is the assay working, averaged over the whole peak set rather than shown
  at one locus. A locus panel can always be cherry-picked; a metaprofile cannot.

PANELS
  A  average profile, forward and reverse, IP and control
  B  the same, split by whether the peak carries the motif -- a motif-positive
     subset with a sharper profile is direct binding; a flat motif-less subset
     is consistent with indirect recruitment
  C  per-peak heatmap sorted by enrichment, so the reader sees the distribution
     rather than only its mean

WHY THE CONTROL IS ON PANEL A
  Without it, any pile-up looks like a footprint. The control is the same
  arithmetic on the no-tag library and should be flat.
"""
import os
import re
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

fwd_bw = dict(zip(snakemake.params.strains, snakemake.input.fwd))     # noqa: F821
rev_bw = dict(zip(snakemake.params.strains, snakemake.input.rev))     # noqa: F821
ctl_bw = dict(zip(snakemake.params.strains, snakemake.input.ctrl))    # noqa: F821
events = dict(zip(snakemake.params.strains, snakemake.input.events))  # noqa: F821
classif = dict(zip(snakemake.params.strains, snakemake.input.classification))  # noqa: F821
repcodes = dict(zip(snakemake.params.strains, snakemake.input.repcodes))  # noqa: F821
STRAINS = list(snakemake.params.strains)                              # noqa: F821
FLANK = int(snakemake.params.flank)                                   # noqa: F821
# Minimum replicate support for a peak to enter the profile. 1 includes every
# called event, which makes this the most general statement about the assay:
# the doublet is visible even across the weakest calls. Raising it to match
# integration.min_reps gives a cleaner profile but a narrower claim.
MIN_REPS = int(snakemake.params.min_reps)                             # noqa: F821
outdir = snakemake.params.outdir                                      # noqa: F821
logp = snakemake.log[0]                                               # noqa: F821

os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

INK, TEAL, RUST = "#1F2933", "#0E6E6B", "#A8471C"
MUTE, GREY, PALE = "#6B7A8A", "#B4BEC7", "#EDF1F3"
CMAP = LinearSegmentedColormap.from_list("perr", ["#FFFFFF", "#C9743F", "#A8471C", "#5C2410"])


def read_points(path):
    """ChExMix event Points, as plain integer coordinates."""
    out = []
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        m = re.match(r"^(?:chr)?(.+):(\d+)$", ln.split("\t")[0])
        if m:
            out.append((m.group(1), int(m.group(2))))
    return out

def read_repcodes(path):
    """peak position -> number of replicates the event was called in."""
    out = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 2 or not f[1].isdigit():
            continue
        m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0])
        if m:
            out[int(m.group(2))] = int(f[1])
    return out

def read_motif_status(path):
    """peak position -> motif-positive, keyed the same way as the events file."""
    st = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            m = re.match(r"^(.+)_(\d+)$", f[0])
            if m:
                st[int(m.group(2))] = (f[i["gbs_perbox_in_window"]] == "True")
    return st


def profile(bw_path, pts, flank):
    """Per-peak signal matrix, rows = peaks, columns = positions.

    Coordinates outside the chromosome are skipped rather than clipped, since a
    clipped row would contribute a truncated profile to the mean.
    """
    bw = pyBigWig.open(bw_path)
    chroms = bw.chroms()
    rows = []
    for chrom, pos in pts:
        c = chrom if chrom in chroms else next(
            (k for k in chroms if k.endswith(chrom) or chrom.endswith(k)), None)
        if c is None:
            continue
        s, e = pos - flank, pos + flank + 1
        if s < 0 or e > chroms[c]:
            continue
        v = np.nan_to_num(np.array(bw.values(c, s, e), dtype=float))
        rows.append(v)
    bw.close()
    return np.vstack(rows) if rows else np.zeros((0, 2 * flank + 1))


x = np.arange(-FLANK, FLANK + 1)
log = open(logp, "w")

fig = plt.figure(figsize=(13, 4.6 * len(STRAINS)))
gs = fig.add_gridspec(len(STRAINS), 3, width_ratios=[1.25, 1.25, 1.0],
                      hspace=0.45, wspace=0.30)

# Colourbars are axes too, so fig.get_axes() would include them in the styling
# loop at the end and undo their own formatting. Track the panel axes instead.
panels = []

for r, S in enumerate(STRAINS):
    pts = read_points(events[S])
    mot = read_motif_status(classif[S])
    if MIN_REPS > 1:
        codes = read_repcodes(repcodes[S])
        n_before = len(pts)
        pts = [p for p in pts if codes.get(p[1], 0) >= MIN_REPS]
        log.write(f"{S}: reproducibility filter min_reps={MIN_REPS}, "
                  f"{n_before} -> {len(pts)} peaks\n")
    F = profile(fwd_bw[S], pts, FLANK)
    R = profile(rev_bw[S], pts, FLANK)
    C = profile(ctl_bw[S], pts, FLANK)
    log.write(f"{S}: {len(pts)} events, {F.shape[0]} usable windows, flank {FLANK} bp\n")

    # ---- A: average profile, IP forward / reverse / control ----------------
    ax = fig.add_subplot(gs[r, 0])
    panels.append(ax)
    ax.axvline(0, color=INK, lw=.6, zorder=1)
    # A dashed line, not a filled area. Filled, a flat control reads as though
    # it dominates the IP, when the IP's genome-wide background is genuinely
    # lower: ~27% of IP reads sit in peaks and are therefore absent from the
    # background, while the control's reads are spread uniformly.
    ax.plot(x, C.mean(0), color=GREY, lw=1.1, ls="--", label="no-tag control", zorder=2)
    ax.plot(x, F.mean(0), color="#0E5A8A", lw=1.6, label="IP forward 5\u2032", zorder=3)
    ax.plot(x, R.mean(0), color=RUST, lw=1.6, label="IP reverse 5\u2032", zorder=3)
    fmax, rmax = x[F.mean(0).argmax()], x[R.mean(0).argmax()]
    ax.annotate("", xy=(fmax, ax.get_ylim()[1] * .92), xytext=(rmax, ax.get_ylim()[1] * .92),
                arrowprops=dict(arrowstyle="<->", color=MUTE, lw=.9))
    ax.text((fmax + rmax) / 2, ax.get_ylim()[1] * .95, f"{abs(rmax - fmax)} bp",
            ha="center", va="bottom", fontsize=8.5, color=MUTE)
    ax.set_title(f"{S}   average footprint", fontsize=11, color=INK,
                 fontweight="bold", loc="left")
    ax.set_xlabel("distance from crosslink (bp)", fontsize=9.5, color=INK)
    ax.set_ylabel("mean 5\u2032 end signal (CPM)", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    log.write(f"   forward max at {fmax:+d} bp, reverse max at {rmax:+d} bp, "
              f"separation {abs(rmax - fmax)} bp\n")

    # ---- B: split by motif status -----------------------------------------
    ax = fig.add_subplot(gs[r, 1])
    panels.append(ax)
    ax.axvline(0, color=INK, lw=.6)
    idx_pos = [i for i, (_, p) in enumerate(pts) if mot.get(p, False)]
    idx_neg = [i for i, (_, p) in enumerate(pts) if not mot.get(p, False)]
    for idx, col, lab in ((idx_pos, RUST, "motif-positive"),
                          (idx_neg, TEAL, "motif-less")):
        if not idx or max(idx) >= F.shape[0]:
            continue
        tot = (F[idx] + R[idx]).mean(0)
        ax.plot(x, tot, color=col, lw=1.7, label=f"{lab}  (n={len(idx)})")
    ax.set_title(f"{S}   motif-positive vs motif-less", fontsize=11, color=INK,
                 fontweight="bold", loc="left")
    ax.set_xlabel("distance from crosslink (bp)", fontsize=9.5, color=INK)
    ax.set_ylabel("mean 5\u2032 end signal, both strands", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    # ---- C: per-peak heatmap ----------------------------------------------
    ax = fig.add_subplot(gs[r, 2])
    panels.append(ax)
    M = F + R
    if M.shape[0]:
        order = np.argsort(-M.sum(1))
        vmax = np.percentile(M, 99) or 1
        im = ax.imshow(M[order], aspect="auto", cmap=CMAP, vmin=0, vmax=vmax,
                       extent=[-FLANK, FLANK, M.shape[0], 0], interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
        cb.set_label("5\u2032 end signal (CPM)", fontsize=8, color=INK)
        cb.ax.tick_params(labelsize=7.5, colors=INK)
        cb.outline.set_visible(False)
    ax.axvline(0, color=INK, lw=.6)
    ax.set_title(f"{S}   every peak, sorted by signal", fontsize=11, color=INK,
                 fontweight="bold", loc="left")
    ax.set_xlabel("distance from crosslink (bp)", fontsize=9.5, color=INK)
    ax.set_ylabel("binding events", fontsize=9.5, color=INK)

for ax in panels:
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GREY)
    ax.tick_params(colors=INK, labelsize=8.5)

fig.suptitle("ChIP-exo footprint architecture, averaged over all binding events",
             fontsize=15, fontweight="bold", color=INK, x=0.02, ha="left", y=0.995)
fig.text(0.02, 0.965,
         "Forward-strand 5\u2032 ends mark the left border of the footprint and "
         "reverse-strand ends the right; their separation is the protected width.",
         fontsize=10, color=MUTE, ha="left")

fig.savefig(os.path.join(outdir, "metaprofile.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(outdir, "metaprofile.png"), dpi=300, bbox_inches="tight")
log.write("wrote metaprofile.pdf\n")
log.close()