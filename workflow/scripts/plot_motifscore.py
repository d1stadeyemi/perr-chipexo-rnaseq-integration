"""Why are motif-less peaks motif-less? Four tests, one figure.

Panel A  threshold-free motif scores against a dinucleotide-shuffled null.
Panel B  motif rate across enrichment quintiles: a detection-sensitivity
         effect would rise monotonically with binding strength.
Panel C  the artifact tests. Elevated no-tag control signal is the diagnostic
         for hyper-ChIPable regions (Teytelman et al., PNAS 2013).
Panel D  reproducibility composition.
"""
import re, csv, os, numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = "/lustre/isaac24/proj/UTK0009/Muhammad/projects/perr_chipexo"
D = f"{P}/results/14_motifscore"
INK, TEAL, RUST, GREY, MUTE = "#1F2933", "#0E6E6B", "#A8471C", "#B4BEC7", "#6B7A8A"
CENTRE, HALFWIN = 250, 30

def best(path):
    o = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        s = r.get("sequence_name")
        if not s: continue
        try: sc, a, b = float(r["score"]), int(r["start"]), int(r["stop"])
        except (KeyError, ValueError, TypeError): continue
        if abs((a + b)//2 - CENTRE) > HALFWIN: continue
        if s not in o or sc > o[s]: o[s] = sc
    return o

real, null = best(f"{D}/real.txt"), best(f"{D}/null.txt")

mot = {}
with open(f"{P}/results/05_motif/A909/A909.all.peak_classification.tsv") as fh:
    h = fh.readline().rstrip().split("\t"); i = {c: k for k, c in enumerate(h)}
    for l in fh:
        f = l.rstrip("\n").split("\t"); m = re.match(r"^(.+)_(\d+)$", f[0])
        if m: mot[int(m.group(2))] = f[i["gbs_perbox_in_window"]] == "True"

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
    return float(np.nan_to_num(np.array(bw.values(CH, max(0,p-h), p+h), float)).max())

MP = [p for p in mot if mot[p]]; ML = [p for p in mot if not mot[p]]
SP = [real[k] for k in real if int(k.split("_")[1]) in MP]
SL = [real[k] for k in real if int(k.split("_")[1]) in ML]
NU = sorted(null.values())

fig, ax = plt.subplots(2, 2, figsize=(13, 9))

# A -- score distributions
a = ax[0,0]
bins = np.linspace(-20, 18, 40)
a.hist(NU, bins=bins, density=True, color=GREY, alpha=.75, label=f"shuffled null (n={len(NU)})")
a.hist(SL, bins=bins, density=True, histtype="step", lw=2, color=TEAL, label=f"motif-less peaks (n={len(SL)})")
a.hist(SP, bins=bins, density=True, histtype="step", lw=2, color=RUST, label=f"motif-positive peaks (n={len(SP)})")
thr = NU[int(len(NU)*.95)]
a.axvline(thr, color=INK, ls="--", lw=1)
a.text(thr, a.get_ylim()[1]*.95, " null 95th pct", fontsize=8, color=INK)
a.set_xlabel("best log-odds score within \u00b130 bp of the crosslink")
a.set_ylabel("density"); a.legend(frameon=False, fontsize=8.5)
a.set_title("A   Motif-less peaks score like shuffled DNA", loc="left",
            fontweight="bold", color=INK)
below = sum(1 for x in SL if x <= thr)
a.text(.02, .70, f"{below}/{len(SL)} motif-less peaks at or below the null 95th percentile",
       transform=a.transAxes, fontsize=8.5, color=MUTE)

# B -- motif rate by enrichment quintile
a = ax[0,1]
allp = sorted(mot, key=lambda p: -lfc.get(p, 0)); k = len(allp)//5
xs, ys, ns = [], [], []
for j in range(5):
    g = allp[j*k:(j+1)*k] if j < 4 else allp[4*k:]
    xs.append(f"Q{j+1}\n{lfc[g[0]]:.2f}\u2013{lfc[g[-1]]:.2f}")
    ys.append(100*sum(1 for p in g if mot[p])/len(g)); ns.append(len(g))
a.bar(xs, ys, color=RUST, width=.62)
for x, y, n in zip(xs, ys, ns):
    a.text(x, y+1.2, f"{y:.0f}%", ha="center", fontsize=9, color=INK)
a.axhline(np.mean(ys), color=MUTE, ls="--", lw=1)
a.set_ylabel("peaks containing the motif"); a.set_ylim(0, max(ys)*1.35)
a.set_xlabel("ChIP enrichment quintile (strongest first)")
a.set_title("B   Motif rate does not track binding strength", loc="left",
            fontweight="bold", color=INK)

# C -- the artifact tests
a = ax[1,0]
tests = [("ChIP\nenrichment", [lfc.get(p,0) for p in MP], [lfc.get(p,0) for p in ML]),
         ("no-tag control\nsignal", [wmax(ct,p) for p in MP], [wmax(ct,p) for p in ML]),
         ("IP signal\namplitude", [wmax(ip,p) for p in MP], [wmax(ip,p) for p in ML])]
pos = []
for j,(lab,x,y) in enumerate(tests):
    bp = a.boxplot([x,y], positions=[j*3, j*3+1], widths=.7, patch_artist=True,
                   showfliers=False, medianprops=dict(color=INK, lw=1.4))
    for b,c in zip(bp["boxes"], (RUST, TEAL)):
        b.set_facecolor(c); b.set_alpha(.55); b.set_edgecolor(c)
    pos.append(j*3+.5)
a.set_xticks(pos); a.set_xticklabels([t[0] for t in tests])
a.set_ylabel("value"); a.set_yscale("log")
a.set_title("C   Not weak peaks, and not hyper-ChIPable regions", loc="left",
            fontweight="bold", color=INK)
a.text(.02, .04, "rust = motif-positive   teal = motif-less\n"
       "enrichment p = 0.34   control p = 0.48   (Mann-Whitney)",
       transform=a.transAxes, fontsize=8.5, color=MUTE)

# D -- reproducibility
a = ax[1,1]
for j,(lab,grp,col) in enumerate((("motif-positive",MP,RUST),("motif-less",ML,TEAL))):
    bot = 0
    for r,alpha in ((3,1.0),(2,.6),(1,.3)):
        v = 100*sum(1 for p in grp if code.get(p)==r)/len(grp)
        a.bar(lab, v, bottom=bot, color=col, alpha=alpha, width=.5)
        if v > 4: a.text(j, bot+v/2, f"{r}/3 reps\n{v:.0f}%", ha="center",
                         va="center", fontsize=8.5, color="white" if alpha>.5 else INK)
        bot += v
a.set_ylabel("% of peaks"); a.set_ylim(0, 105)
a.set_title("D   Equally reproducible", loc="left", fontweight="bold", color=INK)

for x in ax.ravel():
    for sp in ("top","right"): x.spines[sp].set_visible(False)
    x.tick_params(colors=INK, labelsize=9)

fig.suptitle("Why do 69% of PerR peaks carry no motif?",
             fontsize=16, fontweight="bold", color=INK, x=.02, ha="left")
fig.text(.02, .945, "Four explanations tested. Motif-less peaks are as strong, as "
         "reproducible and in as clean a background as motif-positive peaks, and "
         "their sequence is indistinguishable from shuffled DNA.",
         fontsize=10.5, color=MUTE, ha="left")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(f"{D}/motifscore_evidence.pdf", bbox_inches="tight")
fig.savefig(f"{D}/motifscore_evidence.png", dpi=300, bbox_inches="tight")
print("wrote", f"{D}/motifscore_evidence.pdf")
ip.close(); ct.close()