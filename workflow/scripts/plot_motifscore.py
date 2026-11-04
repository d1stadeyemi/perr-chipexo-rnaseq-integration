"""Why are motif-less peaks motif-less? Four tests, four separate figures (A909).

  motifscore_scores           threshold-free motif scores vs a dinucleotide-shuffled null
  motifscore_quintiles        motif rate across ChIP-enrichment quintiles
  motifscore_artifacts        ChIP enrichment, control signal and IP amplitude,
                              with vs without motif (one sub-axis per measure)
  motifscore_reproducibility  replicate support, with vs without motif

No titles or explanatory text. Counts, thresholds and test statistics are
printed to the terminal for the figure legend.
Run:  python workflow/scripts/plot_motifscore.py
"""
import re, csv, numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = "/lustre/isaac24/proj/UTK0009/Muhammad/projects/perr_chipexo"
D = f"{P}/results/14_motifscore"
CENTRE, HALFWIN = 250, 30

# ============================================================================
# FONT AND SIZE CONTROLS -- change these to resize everything
# ============================================================================
TICK_SIZE   = 20     # axis tick numbers and category labels
LABEL_SIZE  = 22     # axis titles
LEGEND_SIZE = 19     # legend text
BAR_TEXT    = 18     # numbers printed on / above bars
FIG_W, FIG_H = 8, 6  # size of each saved figure (inches)
SHOW_P = True        # print Mann-Whitney p above the boxes in the artifact figure
# ============================================================================

# perR locus: plasmid copy-number artifact, excluded as in the Methods
EXCLUDE = (410242, 410922)

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.labelsize": LABEL_SIZE, "xtick.labelsize": TICK_SIZE,
                     "ytick.labelsize": TICK_SIZE})
INK, TEAL, RUST, GREY = "#1F2933", "#0E6E6B", "#A8471C", "#B4BEC7"
LEG = {"weight": "bold", "size": LEGEND_SIZE}


def keep(p):
    return not (EXCLUDE[0] <= p <= EXCLUDE[1])


def best(path):
    o = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        s = r.get("sequence_name")
        if not s: continue
        try: sc, a, b = float(r["score"]), int(r["start"]), int(r["stop"])
        except (KeyError, ValueError, TypeError): continue
        if abs((a + b) // 2 - CENTRE) > HALFWIN: continue
        if s not in o or sc > o[s]: o[s] = sc
    return o


def pos_of(name):
    return int(name.split("_")[1])


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK); ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(colors=INK, width=1.4, length=6)
    for t in ax.get_xticklabels() + ax.get_yticklabels(): t.set_fontweight("bold")


def save(fig, stem):
    fig.savefig(f"{D}/{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{D}/{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig); print("wrote", f"{D}/{stem}.pdf")


try:
    from scipy.stats import mannwhitneyu
    def mw(x, y): return mannwhitneyu(x, y, alternative="two-sided").pvalue
except ImportError:
    def mw(x, y): return float("nan")

# ---------------------------------------------------------------- data ------
real = {k: v for k, v in best(f"{D}/real.txt").items() if keep(pos_of(k))}
null = {k: v for k, v in best(f"{D}/null.txt").items() if keep(pos_of(k))}

mot = {}
with open(f"{P}/results/05_motif/A909/A909.all.peak_classification.tsv") as fh:
    h = fh.readline().rstrip().split("\t"); i = {c: k for k, c in enumerate(h)}
    for l in fh:
        f = l.rstrip("\n").split("\t"); m = re.match(r"^(.+)_(\d+)$", f[0])
        if m and keep(int(m.group(2))):
            mot[int(m.group(2))] = f[i["gbs_perbox_in_window"]] == "True"

lfc, code = {}, {}
for l in open(f"{P}/results/03_peaks/chexmix/A909/A909.events"):
    if l.startswith("#") or not l.strip(): continue
    f = l.split("\t"); m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0])
    if m: lfc[int(m.group(2))] = float(f[3])
for l in open(f"{P}/results/03_peaks/chexmix/A909/A909.replicationcodes.table"):
    f = l.rstrip().split("\t")
    if len(f) < 2 or not f[1].isdigit(): continue
    m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0])
    if m: code[int(m.group(2))] = int(f[1])

BW = f"{P}/results/10_igv/bigwig/full_depth/A909"
ip, ct = pyBigWig.open(f"{BW}/A909.ip.fwd.cpm.bw"), pyBigWig.open(f"{BW}/A909.ctrl.cpm.bw")
CH = [k for k in ct.chroms() if "007432" in k][0]
def wmax(bw, p, h=100):
    return float(np.nan_to_num(np.array(bw.values(CH, max(0, p - h), p + h), float)).max())

MP = [p for p in mot if mot[p]]; ML = [p for p in mot if not mot[p]]
SP = [real[k] for k in real if pos_of(k) in MP]
SL = [real[k] for k in real if pos_of(k) in ML]
NU = sorted(null.values())
print(f"events: with motif {len(MP)}, without motif {len(ML)} (perR locus excluded)")

# ---------------------------------------------------- 1. score distributions
fig, a = plt.subplots(figsize=(FIG_W, FIG_H))
bins = np.linspace(-20, 18, 40)
a.hist(NU, bins=bins, density=True, color=GREY, alpha=.8, label="Shuffled")
a.hist(SL, bins=bins, density=True, histtype="step", lw=3, color=TEAL, label="Without motif")
a.hist(SP, bins=bins, density=True, histtype="step", lw=3, color=RUST, label="With motif")
thr = NU[int(len(NU) * .95)]
a.axvline(thr, color=INK, ls="--", lw=1.6)
a.set_xlabel("Best motif log-odds score"); a.set_ylabel("Density")
a.legend(frameon=False, prop=LEG, loc="upper left")
style(a); save(fig, "motifscore_scores")
print(f"scores: null 95th percentile = {thr:.2f}; "
      f"{sum(1 for x in SL if x <= thr)}/{len(SL)} events without motif at or below it; "
      f"medians with={np.median(SP):.2f} without={np.median(SL):.2f} null={np.median(NU):.2f}")

# ------------------------------------------- 2. motif rate by enrichment quintile
fig, a = plt.subplots(figsize=(FIG_W, FIG_H))
allp = sorted(mot, key=lambda p: -lfc.get(p, 0)); k = len(allp) // 5
xs, ys = [], []
for j in range(5):
    g = allp[j * k:(j + 1) * k] if j < 4 else allp[4 * k:]
    xs.append(f"Q{j+1}"); ys.append(100 * sum(1 for p in g if mot[p]) / len(g))
    print(f"quintile Q{j+1}: n={len(g)}, log2FC {lfc[g[0]]:.2f} to {lfc[g[-1]]:.2f}, "
          f"motif {ys[-1]:.0f}%")
a.bar(xs, ys, color=RUST, width=.62)
for x, y in zip(xs, ys):
    a.text(x, y + 1.2, f"{y:.0f}%", ha="center", fontsize=BAR_TEXT, fontweight="bold", color=INK)
a.axhline(np.mean(ys), color=INK, ls="--", lw=1.4)
a.set_ylim(0, max(ys) * 1.3)
a.set_ylabel("Events with motif (%)")
a.set_xlabel("ChIP enrichment quintile")
style(a); save(fig, "motifscore_quintiles")

# --------------------------------------- 3. artifact tests, one axis per measure
tests = [("ChIP enrichment\n(log₂ fold)", [lfc.get(p, 0) for p in MP], [lfc.get(p, 0) for p in ML]),
         ("Control signal\n(CPM)",          [wmax(ct, p) for p in MP], [wmax(ct, p) for p in ML]),
         ("IP signal\n(CPM)",               [wmax(ip, p) for p in MP], [wmax(ip, p) for p in ML])]
fig, axs = plt.subplots(1, 3, figsize=(FIG_W * 1.6, FIG_H))
for a, (lab, x, y) in zip(axs, tests):
    bp = a.boxplot([x, y], positions=[0, 1], widths=.6, patch_artist=True, showfliers=False,
                   medianprops=dict(color=INK, lw=2.2),
                   whiskerprops=dict(lw=1.6), capprops=dict(lw=1.6))
    for b, c in zip(bp["boxes"], (RUST, TEAL)):
        b.set_facecolor(c); b.set_alpha(.6); b.set_edgecolor(c); b.set_linewidth(1.6)
    a.set_xticks([0, 1]); a.set_xticklabels(["With", "Without"])
    a.set_ylabel(lab)
    p = mw(x, y)
    print(f"artifact test {lab.split(chr(10))[0]}: median with={np.median(x):.2f} "
          f"without={np.median(y):.2f}, Mann-Whitney p={p:.2g}")
    if SHOW_P and p == p:
        a.text(.5, 1.02, f"p = {p:.2f}", transform=a.transAxes, ha="center",
               fontsize=BAR_TEXT, fontweight="bold", color=INK)
    style(a)
fig.supxlabel("Motif", fontsize=LABEL_SIZE, fontweight="bold")
fig.tight_layout()
save(fig, "motifscore_artifacts")

# --------------------------------------------------------- 4. reproducibility
fig, a = plt.subplots(figsize=(FIG_W * .8, FIG_H))
for j, (lab, grp, col) in enumerate((("With motif", MP, RUST), ("Without motif", ML, TEAL))):
    bot = 0
    for r, alpha in ((3, 1.0), (2, .6), (1, .3)):
        v = 100 * sum(1 for p in grp if code.get(p) == r) / len(grp)
        a.bar(lab, v, bottom=bot, color=col, alpha=alpha, width=.55)
        if v > 6:
            a.text(j, bot + v / 2, f"{r}/3", ha="center", va="center",
                   fontsize=BAR_TEXT, fontweight="bold",
                   color="white" if alpha > .5 else INK)
        print(f"reproducibility {lab}: {r}/3 replicates {v:.0f}%")
        bot += v
a.set_ylabel("Events (%)"); a.set_ylim(0, 105)
style(a); save(fig, "motifscore_reproducibility")

ip.close(); ct.close()
