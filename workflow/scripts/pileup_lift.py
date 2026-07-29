"""METHOD (b): whole-genome-alignment coordinate lift.

Lift each strain-A ChExMix event coordinate into strain-B coordinates using a
nucmer 1-to-1 alignment, then measure strain-B ChIP-exo coverage in a window
centred on the lifted base -- again regardless of whether B called a peak.

WHY THIS EXISTS ALONGSIDE METHOD (a)
  Method (a) anchors on the TSS of an orthologous GENE. That is robust but
  blind to any event with no gene assignment, and it cannot address intergenic
  sites precisely -- including the strongest site in this dataset, the 152 bp
  divergent intergenic region at A909 874,861 between queA and the VIT1/CCC1
  transporter, where a bound repressor occludes BOTH promoters. This method
  lifts the exact base and so covers those cases.

  It is also a genuinely INDEPENDENT line of evidence: it uses DNA alignment
  rather than protein reciprocal best hits. Agreement between two methods that
  share no intermediate is worth far more than either alone, which is why the
  two outputs are kept in separate files and never averaged.

THE NULL
  Signal at a lifted event position is meaningless without knowing what an
  arbitrary position looks like. `n_null_positions` random liftable positions
  are sampled with a fixed seed and measured identically. Random genomic
  positions are the correct null here rather than shuffled peaks, because the
  question is whether the orthologous position carries occupancy above the
  genome-wide baseline.

COORDINATE CONVENTIONS -- the easy place to be silently wrong
  * ChExMix Points are 1-based, and prefixed 'chr' in the events file.
  * nucmer show-coords output is 1-based inclusive on both genomes.
  * A reversed alignment block has S2 > E2; the lift must then run backwards.
  * pyBigWig queries are 0-based half-open.
  Each conversion is done explicitly below rather than assumed.
"""
import os
import random
import re
import statistics

import pyBigWig

coords = snakemake.input.coords              # noqa: F821
events_a = snakemake.input.events_a          # noqa: F821
events_b = snakemake.input.events_b          # noqa: F821
bw_a_ip_p = snakemake.input.bw_a_ip          # noqa: F821
bw_b_ip_p = snakemake.input.bw_b_ip          # noqa: F821
bw_b_ct_p = snakemake.input.bw_b_ct          # noqa: F821
out = snakemake.output.tsv                   # noqa: F821
summ = snakemake.output.summary              # noqa: F821
sa = snakemake.params.strain_a               # noqa: F821
sb = snakemake.params.strain_b               # noqa: F821
HALF = int(snakemake.params.halfwin)         # noqa: F821
N_NULL = int(snakemake.params.n_null)        # noqa: F821
SEED = int(snakemake.params.seed)            # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


# ---------------------------------------------------------------------------
# alignment blocks
# ---------------------------------------------------------------------------
def read_coords(path):
    """show-coords -T -H -r  ->  S1 E1 S2 E2 LEN1 LEN2 %IDY [..] TAGR TAGQ

    Returns [(rchrom, s1, e1, qchrom, s2, e2), ...] with s1 <= e1. A reversed
    block keeps s2 > e2, which the lift below detects and handles."""
    blocks = []
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            try:
                s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            except ValueError:
                continue
            rchrom, qchrom = f[-2], f[-1]
            if s1 > e1:                      # reference span is always forward
                s1, e1 = e1, s1
                s2, e2 = e2, s2
            blocks.append((rchrom, s1, e1, qchrom, s2, e2))
    return blocks


BLOCKS = read_coords(coords)


def lift(rchrom, pos):
    """Lift a 1-based reference position into 1-based query coordinates.
    Returns (qchrom, qpos, orientation) or None if the position is not covered
    by any alignment block -- which is itself informative and is reported."""
    for (rc, s1, e1, qc, s2, e2) in BLOCKS:
        if rc != rchrom or not (s1 <= pos <= e1):
            continue
        off = pos - s1
        if s2 <= e2:
            return (qc, s2 + off, "+")
        return (qc, s2 - off, "-")
    return None


# ---------------------------------------------------------------------------
# events
# ---------------------------------------------------------------------------
def read_events(path):
    """-> [(chrom_without_prefix, 1-based position), ...]"""
    ev = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            point = ln.split("\t")[0].strip()
            m = re.match(r"^(?:chr)?(.+):(\d+)$", point)
            if m:
                ev.append((m.group(1), int(m.group(2))))
    return ev


EV_A, EV_B = read_events(events_a), read_events(events_b)


# ---------------------------------------------------------------------------
# bigwigs
# ---------------------------------------------------------------------------
def open_bw(path):
    bw = pyBigWig.open(path)
    if bw is None:
        raise SystemExit(f"could not open bigwig: {path}")
    return bw


BWA, BWB, BWBC = open_bw(bw_a_ip_p), open_bw(bw_b_ip_p), open_bw(bw_b_ct_p)


def resolve_chrom(bw, name):
    chroms = bw.chroms()
    for cand in (name, "chr" + name, name.replace("chr", "", 1)):
        if cand in chroms:
            return cand, chroms[cand]
    raise SystemExit(
        f"sequence '{name}' not found in bigwig. Available: {list(chroms)[:5]}")


def stats_at(bw, chrom, pos1, half):
    """Mean/max coverage in +/- half bp around a 1-based position."""
    c, clen = resolve_chrom(bw, chrom)
    s = max(0, pos1 - 1 - half)
    e = min(clen, pos1 - 1 + half + 1)
    if e <= s:
        return ("", "")
    try:
        mean = bw.stats(c, s, e, type="mean")[0]
        mx = bw.stats(c, s, e, type="max")[0]
    except (RuntimeError, IndexError):
        return ("", "")
    return (0.0 if mean is None else mean, 0.0 if mx is None else mx)


# strain-B events, for "did B independently call a peak near the lifted base?"
B_POS = {}
for c, p in EV_B:
    B_POS.setdefault(c, []).append(p)
for c in B_POS:
    B_POS[c].sort()


def b_called_within(qchrom, qpos, tol):
    for p in B_POS.get(qchrom, []):
        if abs(p - qpos) <= tol:
            return p
    return ""


# ---------------------------------------------------------------------------
# lift every A event
# ---------------------------------------------------------------------------
rows, sig_event, n_unlifted = [], [], 0

for chrom, pos in sorted(EV_A):
    a_mean, a_max = stats_at(BWA, chrom, pos, HALF)
    L = lift(chrom, pos)
    if L is None:
        n_unlifted += 1
        rows.append([f"{chrom}:{pos}", chrom, pos, a_mean, a_max,
                     "", "", "", "", "", "", "", "UNLIFTABLE"])
        continue
    qc, qp, orient = L
    b_mean, b_max = stats_at(BWB, qc, qp, HALF)
    c_mean, c_max = stats_at(BWBC, qc, qp, HALF)
    near = b_called_within(qc, qp, HALF * 2)
    ratio = ""
    try:
        if float(c_max) > 0:
            ratio = round(float(b_max) / float(c_max), 4)
    except (TypeError, ValueError):
        pass
    rows.append([f"{chrom}:{pos}", chrom, pos, a_mean, a_max,
                 qc, qp, orient, b_mean, b_max, c_mean, ratio,
                 "B_PEAK_NEARBY" if near != "" else "B_NO_PEAK"])
    if b_max != "":
        sig_event.append(float(b_max))

# ---------------------------------------------------------------------------
# matched random null
# ---------------------------------------------------------------------------
rng = random.Random(SEED)
ref_chrom = EV_A[0][0] if EV_A else None
ref_len = None
if ref_chrom:
    _, ref_len = resolve_chrom(BWA, ref_chrom)

null_vals, n_null_done, attempts = [], 0, 0
while ref_len and n_null_done < N_NULL and attempts < N_NULL * 20:
    attempts += 1
    p = rng.randint(1, ref_len)
    L = lift(ref_chrom, p)
    if L is None:
        continue
    qc, qp, _ = L
    _, b_max = stats_at(BWB, qc, qp, HALF)
    if b_max == "":
        continue
    null_vals.append(float(b_max))
    n_null_done += 1

HDR = ["event", f"{sa}_chrom", f"{sa}_pos", f"{sa}_ip_mean", f"{sa}_ip_max",
       f"{sb}_chrom", f"{sb}_pos", "orientation",
       f"{sb}_ip_mean", f"{sb}_ip_max", f"{sb}_ctrl_mean",
       f"{sb}_ip_over_ctrl", "status"]

with open(out, "w") as ob:
    ob.write("\t".join(HDR) + "\n")
    for r in rows:
        ob.write("\t".join(str(v) for v in r) + "\n")


def describe(vals):
    v = sorted(vals)
    if not v:
        return dict(n=0, median="", q75="", q90="", mean="")
    n = len(v)

    def q(p):
        return v[min(n - 1, int(p * n))]

    return dict(n=n, median=round(statistics.median(v), 4),
                q75=round(q(0.75), 4), q90=round(q(0.90), 4),
                mean=round(statistics.fmean(v), 4))


d_ev, d_nl = describe(sig_event), describe(null_vals)

with open(summ, "w") as sf:
    sf.write("group\tn\tmedian_B_ip_max\tq75\tq90\tmean\n")
    for nm, d in (("LIFTED_A_EVENTS", d_ev), ("RANDOM_NULL", d_nl)):
        sf.write(f"{nm}\t{d['n']}\t{d['median']}\t{d['q75']}\t{d['q90']}\t{d['mean']}\n")

with open(logp, "w") as lg:
    lg.write(f"alignment blocks: {len(BLOCKS)}\n")
    lg.write(f"{sa} events: {len(EV_A)}   unliftable: {n_unlifted}\n")
    lg.write(f"window: +/- {HALF} bp around the lifted base\n")
    lg.write(f"null positions sampled: {d_nl['n']} (seed {SEED})\n")
    lg.write(f"lifted events   : median={d_ev['median']} q90={d_ev['q90']}\n")
    lg.write(f"random null     : median={d_nl['median']} q90={d_nl['q90']}\n")
    try:
        lg.write(f"\nevent median / null median = "
                 f"{float(d_ev['median']) / float(d_nl['median']):.2f}\n")
        lg.write("  >> clearly above 1: orthologous positions carry occupancy "
                 f"in {sb} even where no peak was called -> SENSITIVITY\n")
        lg.write("  >> near 1: those positions are genuinely unbound -> BIOLOGY\n")
    except (ZeroDivisionError, ValueError, TypeError):
        pass

for bw in (BWA, BWB, BWBC):
    bw.close()