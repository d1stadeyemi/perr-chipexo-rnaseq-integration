"""METHOD (a): ortholog-anchored signal extraction.

For every gene bound by PerR in strain A, find its reciprocal-best-hit ortholog
in strain B and measure the ChIP-exo coverage in a window at that orthologous
gene's TSS -- REGARDLESS of whether strain B called a peak there.

The comparison that matters is against a NULL. Signal is therefore also
extracted at the orthologous TSS of every 1:1 ortholog whose strain-A gene is
NOT bound. That distribution defines what "background" means in these units,
so "elevated but sub-threshold" becomes a quantitative statement rather than an
impression.

Groups reported in `group`:
  BOUND_BOTH   bound in A and in B          -> positive control; these SHOULD
                                               show clear signal in B
  A_ONLY       bound in A, not called in B  -> THE SET UNDER TEST
  B_ONLY       bound in B, not called in A  -> reverse direction
  NULL         bound in neither             -> background distribution

READING THE RESULT
  If A_ONLY sites in strain B sit clearly above NULL, the sites are occupied
  and the peak caller lacked power -- the half-slurry sensitivity explanation.
  If A_ONLY is indistinguishable from NULL, the sites are genuinely unbound in
  strain B -- real cistrome divergence.

CAVEAT built into this method: it anchors on the TSS of the ORTHOLOGOUS GENE,
not on the exact orthologous base. Promoter architecture can shift between
strains, so the window must be wide enough to tolerate that (see
config pileup.tss_upstream / tss_downstream). Events with no gene assignment --
including intergenic sites -- are invisible to this method entirely; that is
what method (b) is for.
"""
import os
import re
import statistics

import pyBigWig

orthologs = snakemake.input.orthologs        # noqa: F821
assign_a = snakemake.input.assign_a          # noqa: F821
assign_b = snakemake.input.assign_b          # noqa: F821
gff_a = snakemake.input.gff_a                # noqa: F821
gff_b = snakemake.input.gff_b                # noqa: F821
bw_a_ip = snakemake.input.bw_a_ip            # noqa: F821
bw_a_ct = snakemake.input.bw_a_ct            # noqa: F821
bw_b_ip = snakemake.input.bw_b_ip            # noqa: F821
bw_b_ct = snakemake.input.bw_b_ct            # noqa: F821
out = snakemake.output.tsv                   # noqa: F821
summ = snakemake.output.summary              # noqa: F821
sa = snakemake.params.strain_a               # noqa: F821
sb = snakemake.params.strain_b               # noqa: F821
UP = int(snakemake.params.upstream)          # noqa: F821
DOWN = int(snakemake.params.downstream)      # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


# ---------------------------------------------------------------------------
# genes
# ---------------------------------------------------------------------------
def gene_coords(gff):
    """locus_tag -> (chrom, start, end, strand). 1-based inclusive, as in GFF."""
    g = {}
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            lt = re.search(r"locus_tag=([^;]+)", f[8])
            if lt:
                g[lt.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return g


def tss(rec):
    """Transcription start = gene start on +, gene end on -."""
    _, s, e, strand = rec
    return s if strand == "+" else e


GA, GB = gene_coords(gff_a), gene_coords(gff_b)


def bound_genes(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = hdr.index("gene_id")
        return {ln.rstrip("\n").split("\t")[i] for ln in fh if ln.strip()}


BA, BB = bound_genes(assign_a), bound_genes(assign_b)

# 1:1 orthologs only -- ambiguous pairs cannot be anchored unambiguously
pairs = []
with open(orthologs) as fh:
    fh.readline()
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 7 and f[6] == "1:1":
            pairs.append((f[0], f[1]))


# ---------------------------------------------------------------------------
# bigwig access
# ---------------------------------------------------------------------------
def open_bw(path):
    bw = pyBigWig.open(path)
    if bw is None:
        raise SystemExit(f"could not open bigwig: {path}")
    return bw


BW = {k: open_bw(p) for k, p in
      (("a_ip", bw_a_ip), ("a_ct", bw_a_ct), ("b_ip", bw_b_ip), ("b_ct", bw_b_ct))}


def resolve_chrom(bw, name):
    """BAM/bigwig sequence names may or may not carry a 'chr' prefix relative to
    the GFF. Resolve rather than assume, and fail loudly if neither form is
    present -- a silent mismatch here would return zeros everywhere and look
    exactly like 'no signal'."""
    chroms = bw.chroms()
    for cand in (name, "chr" + name, name.replace("chr", "", 1)):
        if cand in chroms:
            return cand, chroms[cand]
    raise SystemExit(
        f"sequence '{name}' not found in bigwig. Available: {list(chroms)[:5]}")


def window_stats(bw, chrom, centre, up, down, strand):
    """Mean and max coverage in a strand-aware window around a TSS.
    Returns (mean, max, start, end) in 0-based half-open coordinates."""
    c, clen = resolve_chrom(bw, chrom)
    if strand == "+":
        s, e = centre - up, centre + down
    else:
        s, e = centre - down, centre + up
    s = max(0, s)
    e = min(clen, e)
    if e <= s:
        return ("", "", s, e)
    try:
        mean = bw.stats(c, s, e, type="mean")[0]
        mx = bw.stats(c, s, e, type="max")[0]
    except (RuntimeError, IndexError):
        return ("", "", s, e)
    return (0.0 if mean is None else mean, 0.0 if mx is None else mx, s, e)


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------
rows = []
groups = {}

for ga, gb in pairs:
    if ga not in GA or gb not in GB:
        continue
    in_a, in_b = ga in BA, gb in BB
    if in_a and in_b:
        grp = "BOUND_BOTH"
    elif in_a:
        grp = "A_ONLY"
    elif in_b:
        grp = "B_ONLY"
    else:
        grp = "NULL"

    ca, sa_, ea_, stra = GA[ga]
    cb, sb_, eb_, strb = GB[gb]
    ta, tb = tss(GA[ga]), tss(GB[gb])

    a_ip_mean, a_ip_max, _, _ = window_stats(BW["a_ip"], ca, ta, UP, DOWN, stra)
    a_ct_mean, a_ct_max, _, _ = window_stats(BW["a_ct"], ca, ta, UP, DOWN, stra)
    b_ip_mean, b_ip_max, bs, be = window_stats(BW["b_ip"], cb, tb, UP, DOWN, strb)
    b_ct_mean, b_ct_max, _, _ = window_stats(BW["b_ct"], cb, tb, UP, DOWN, strb)

    def ratio(ip, ct):
        try:
            return round(float(ip) / float(ct), 4) if float(ct) > 0 else ""
        except (TypeError, ValueError):
            return ""

    rows.append([ga, gb, grp,
                 cb, bs, be, strb,
                 a_ip_mean, a_ip_max, a_ct_mean,
                 b_ip_mean, b_ip_max, b_ct_mean,
                 ratio(b_ip_max, b_ct_max)])
    groups.setdefault(grp, []).append(b_ip_max if b_ip_max != "" else 0.0)

HDR = [f"{sa}_gene", f"{sb}_gene", "group",
       f"{sb}_chrom", f"{sb}_win_start", f"{sb}_win_end", f"{sb}_strand",
       f"{sa}_ip_mean", f"{sa}_ip_max", f"{sa}_ctrl_mean",
       f"{sb}_ip_mean", f"{sb}_ip_max", f"{sb}_ctrl_mean",
       f"{sb}_ip_over_ctrl"]

with open(out, "w") as ob:
    ob.write("\t".join(HDR) + "\n")
    for r in sorted(rows, key=lambda x: (x[2], x[0])):
        ob.write("\t".join(str(v) for v in r) + "\n")


# ---------------------------------------------------------------------------
# summarise: is A_ONLY above NULL in strain B?
# ---------------------------------------------------------------------------
def describe(vals):
    v = sorted(float(x) for x in vals if x != "")
    if not v:
        return dict(n=0, median="", q75="", q90="", mean="")
    n = len(v)

    def q(p):
        return v[min(n - 1, int(p * n))]

    return dict(n=n, median=round(statistics.median(v), 4),
                q75=round(q(0.75), 4), q90=round(q(0.90), 4),
                mean=round(statistics.fmean(v), 4))


with open(summ, "w") as sf:
    sf.write("group\tn\tmedian_B_ip_max\tq75\tq90\tmean\n")
    for g in ("BOUND_BOTH", "A_ONLY", "B_ONLY", "NULL"):
        d = describe(groups.get(g, []))
        sf.write(f"{g}\t{d['n']}\t{d['median']}\t{d['q75']}\t{d['q90']}\t{d['mean']}\n")

with open(logp, "w") as lg:
    lg.write(f"1:1 ortholog pairs considered: {len(pairs)}\n")
    lg.write(f"rows written: {len(rows)}\n")
    lg.write(f"window: TSS -{UP} / +{DOWN} bp (strand-aware)\n")
    for g in ("BOUND_BOTH", "A_ONLY", "B_ONLY", "NULL"):
        d = describe(groups.get(g, []))
        lg.write(f"  {g}: n={d['n']} median={d['median']} q90={d['q90']}\n")
    nul = describe(groups.get("NULL", []))
    aon = describe(groups.get("A_ONLY", []))
    if nul["n"] and aon["n"] and nul["median"]:
        try:
            lg.write(f"\nA_ONLY median / NULL median = "
                     f"{float(aon['median']) / float(nul['median']):.2f}\n")
            lg.write("  >> clearly above 1 supports SENSITIVITY (sites occupied, "
                     "peak caller underpowered)\n")
            lg.write("  >> near 1 supports BIOLOGY (sites genuinely unbound in "
                     f"{sb})\n")
        except (ZeroDivisionError, ValueError):
            pass

for bw in BW.values():
    bw.close()