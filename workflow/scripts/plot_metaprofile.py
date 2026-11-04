"""plot_metaprofile.py -- the canonical ChIP-exo footprint figures.

Writes THREE SEPARATE figures per strain (no titles, no explanatory text):

  metaprofile_footprint_<strain>   average forward / reverse 5' end signal and
                                   the untagged control, aligned on every event
  metaprofile_motif_<strain>       the same signal (both strands summed), split by
                                   whether the event carries the GBS motif
  metaprofile_heatmap_<strain>     per-event heatmap sorted by signal

Each is written as .pdf (vector, for the paper) and .png (300 dpi, for slides).
Panel letters, strain labels and captions are added in the figure layout.
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
MIN_REPS = int(snakemake.params.min_reps)                             # noqa: F821
outdir = snakemake.params.outdir                                      # noqa: F821
logp = snakemake.log[0]                                               # noqa: F821

os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

# ============================================================================
# FONT AND SIZE CONTROLS -- change these to resize everything
# ============================================================================
TICK_SIZE   = 20     # axis tick numbers
LABEL_SIZE  = 22     # axis titles
LEGEND_SIZE = 19     # legend text
ANNOT_SIZE  = 20     # the "44 bp" footprint label
CBAR_SIZE   = 18     # heatmap colour-bar label and ticks
LINE_W      = 3.0    # width of the profile lines
FIG_W, FIG_H = 8, 6  # size of each saved figure (inches)
# ============================================================================

# Events excluded from the profiles: the perR locus, where the expression
# plasmid inflates read depth independent of binding (see Methods).
# Add the CJB111 perR window here if its events fall there.
EXCLUDE = {"A909": ("NC_007432.1", 410242, 410922)}

plt.rcParams.update({
    "font.weight": "bold", "axes.labelweight": "bold",
    "axes.labelsize": LABEL_SIZE, "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE, "legend.fontsize": LEGEND_SIZE,
})

INK, TEAL, RUST = "#1F2933", "#0E6E6B", "#A8471C"
GREY, BLUE = "#B4BEC7", "#0E5A8A"
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
    Windows running off the chromosome are skipped rather than clipped.
    Also returns the points actually used, so row indices always line up
    with the motif labels."""
    bw = pyBigWig.open(bw_path)
    chroms = bw.chroms()
    rows, used = [], []
    for chrom, pos in pts:
        c = chrom if chrom in chroms else next(
            (k for k in chroms if k.endswith(chrom) or chrom.endswith(k)), None)
        if c is None:
            continue
        s, e = pos - flank, pos + flank + 1
        if s < 0 or e > chroms[c]:
            continue
        rows.append(np.nan_to_num(np.array(bw.values(c, s, e), dtype=float)))
        used.append((chrom, pos))
    bw.close()
    M = np.vstack(rows) if rows else np.zeros((0, 2 * flank + 1))
    return M, used


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK)
        ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(colors=INK, width=1.4, length=6)
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontweight("bold")


def save(fig, stem):
    fig.savefig(os.path.join(outdir, stem + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(outdir, stem + ".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    log.write(f"   wrote {stem}\n")


x = np.arange(-FLANK, FLANK + 1)
XLAB = "Distance from binding event centre (bp)"
LEG = {"weight": "bold", "size": LEGEND_SIZE}
log = open(logp, "w")

for S in STRAINS:
    pts = read_points(events[S])
    if S in EXCLUDE:
        c0, a, b = EXCLUDE[S]
        n0 = len(pts)
        pts = [p for p in pts if not (p[0].endswith(c0) and a <= p[1] <= b)]
        log.write(f"{S}: excluded {n0 - len(pts)} events in {c0}:{a}-{b}\n")
    if MIN_REPS > 1:
        codes = read_repcodes(repcodes[S])
        n0 = len(pts)
        pts = [p for p in pts if codes.get(p[1], 0) >= MIN_REPS]
        log.write(f"{S}: min_reps={MIN_REPS}, {n0} -> {len(pts)} events\n")
    mot = read_motif_status(classif[S])
    F, used = profile(fwd_bw[S], pts, FLANK)
    R, _ = profile(rev_bw[S], pts, FLANK)
    C, _ = profile(ctl_bw[S], pts, FLANK)
    log.write(f"{S}: {F.shape[0]} usable windows, flank {FLANK} bp\n")

    # ---- 1. average footprint ----------------------------------------------
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color=INK, lw=1.0, zorder=1)
    ax.plot(x, C.mean(0), color=GREY, lw=LINE_W * .8, ls="--", label="Control", zorder=2)
    ax.plot(x, F.mean(0), color=BLUE, lw=LINE_W, label="+ strand", zorder=3)
    ax.plot(x, R.mean(0), color=RUST, lw=LINE_W, label="− strand", zorder=3)
    fmax, rmax = x[F.mean(0).argmax()], x[R.mean(0).argmax()]
    top = ax.get_ylim()[1]
    ax.annotate("", xy=(fmax, top * .90), xytext=(rmax, top * .90),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.6))
    ax.text((fmax + rmax) / 2, top * .93, f"{abs(rmax - fmax)} bp",
            ha="center", va="bottom", fontsize=ANNOT_SIZE, fontweight="bold", color=INK,
            bbox=dict(facecolor="white", edgecolor="none", pad=2), zorder=4)
    ax.set_xlabel(XLAB)
    ax.set_ylabel("Mean 5′ end signal (CPM)")
    ax.legend(frameon=False, loc="upper right", prop=LEG)
    style(ax)
    save(fig, f"metaprofile_footprint_{S}")
    log.write(f"   forward max {fmax:+d}, reverse max {rmax:+d}, "
              f"separation {abs(rmax - fmax)} bp\n")

    # ---- 2. with vs without motif (counts go to the log, not the legend) ---
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axvline(0, color=INK, lw=1.0)
    idx_pos = [i for i, (_, p) in enumerate(used) if mot.get(p, False)]
    idx_neg = [i for i, (_, p) in enumerate(used) if not mot.get(p, False)]
    for idx, col, lab in ((idx_pos, RUST, "With motif"), (idx_neg, TEAL, "Without motif")):
        if idx:
            ax.plot(x, (F[idx] + R[idx]).mean(0), color=col, lw=LINE_W, label=lab)
            log.write(f"   {lab}: n={len(idx)}\n")
    ax.set_xlabel(XLAB)
    ax.set_ylabel("Mean 5′ end signal (CPM)")
    ax.legend(frameon=False, loc="upper right", prop=LEG)
    style(ax)
    save(fig, f"metaprofile_motif_{S}")

    # ---- 3. per-event heatmap ----------------------------------------------
    fig, ax = plt.subplots(figsize=(FIG_W * .8, FIG_H))
    M = F + R
    if M.shape[0]:
        order = np.argsort(-M.sum(1))
        vmax = np.percentile(M, 99) or 1
        im = ax.imshow(M[order], aspect="auto", cmap=CMAP, vmin=0, vmax=vmax,
                       extent=[-FLANK, FLANK, M.shape[0], 0], interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=.05, pad=.03)
        cb.set_label("5′ end signal (CPM)", fontsize=CBAR_SIZE, fontweight="bold", color=INK)
        cb.ax.tick_params(labelsize=CBAR_SIZE, colors=INK)
        for t in cb.ax.get_yticklabels():
            t.set_fontweight("bold")
        cb.outline.set_visible(False)
    ax.axvline(0, color=INK, lw=1.0)
    ax.set_xlabel(XLAB)
    ax.set_ylabel("Binding events")
    style(ax)
    save(fig, f"metaprofile_heatmap_{S}")

log.close()
