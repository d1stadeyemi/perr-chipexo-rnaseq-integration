"""sens_summary.py

Collapse one sensitivity axis into a table: one row per parameter combination,
with every headline number recomputed from that combination's own outputs.

WHY THE COUNTS ARE RECOMPUTED HERE RATHER THAN READ FROM regulon_summary.tsv
  regulon_summary.tsv counts (peak, gene) ROWS. Every number quoted in the deck
  is UNIQUE GENES, de-duplicated by class precedence. Reading the row counts
  would silently inflate every figure -- which is exactly the confusion that
  produced the "CJB111 has ~13 direct targets" error (rows) when the answer is
  8 (unique genes). The precedence used here is identical to the one in
  build_regulon_table.py and plot_regulon.R.

AXES SUMMARISED AT THEIR OWN LEVEL
  meme     motif width, E-value, palindromy -- the half-site question
  chexmix  peak counts only
  All others report the full downstream chain.
"""
import os
import re

axis = snakemake.params.axis                 # noqa: F821
combos = list(snakemake.params.combos)       # noqa: F821
grid = dict(snakemake.params.grid)           # noqa: F821
SENS = snakemake.params.sensdir              # noqa: F821
strains = list(snakemake.params.strains)     # noqa: F821
out = snakemake.output.tsv                   # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

PREC = ["DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "BOUND_NO_RNA", "INDIRECT"]


def unique_gene_counts(path):
    """Collapse a regulon table to one row per gene using class precedence,
    then count. This is the number that appears in the manuscript."""
    if not os.path.exists(path):
        return {}
    best = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        if "gene_id" not in idx or "classification" not in idx:
            return {}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                f += [""] * (len(hdr) - len(f))
            g, cls = f[idx["gene_id"]], f[idx["classification"]]
            r = PREC.index(cls) if cls in PREC else 99
            if g not in best or r < best[g][0]:
                best[g] = (r, cls)
    counts = {}
    for _, cls in best.values():
        counts[cls] = counts.get(cls, 0) + 1
    counts["_total_genes"] = len(best)
    return counts


def crossstrain(path):
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3:
                d[f[1]] = f[2]
    return d


def rbh_count(path):
    if not os.path.exists(path):
        return ""
    n = 0
    with open(path) as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 7 and f[6] == "1:1":
                n += 1
    return n


def n_events(path):
    if not os.path.exists(path):
        return ""
    return sum(1 for ln in open(path) if not ln.startswith("#") and ln.strip())


def meme_motifs(path):
    """Parse MEME-ChIP summary.tsv: (n_motifs, best consensus, width, E-value)."""
    if not os.path.exists(path):
        return (0, "", "", "")
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            rows.append(f)
    if not rows:
        return (0, "", "", "")
    top = rows[0]
    return (len(rows),
            top[i.get("CONSENSUS", 4)] if len(top) > i.get("CONSENSUS", 4) else "",
            top[i.get("WIDTH", 5)] if len(top) > i.get("WIDTH", 5) else "",
            top[i.get("E-VALUE", 7)] if len(top) > i.get("E-VALUE", 7) else "")


def revcomp(s):
    return s.translate(str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                                     "TGCAYRMKSWVHDBNtgcayrmkswvhdbn"))[::-1]


def palindromic(seq):
    """Crude test: does the consensus resemble its own reverse complement?
    A perfect inverted repeat would score 1.0. Reported as a hint that a wider
    window has recovered a dimer site, NOT as a formal statistic."""
    s = re.sub(r"[^A-Za-z]", "", seq or "")
    if len(s) < 6:
        return ""
    rc = revcomp(s)
    match = sum(1 for a, b in zip(s.upper(), rc.upper()) if a == b)
    return f"{match/len(s):.2f}"


rows, param_keys = [], []
if grid:
    param_keys = sorted(next(iter(grid.values())).keys())

for combo in combos:
    d = os.path.join(SENS, axis, combo)
    rec = {"axis": axis, "combo": combo}
    for k in param_keys:
        rec[k] = grid.get(combo, {}).get(k, "")

    if axis == "meme":
        for s in strains:
            n, cons, w, e = meme_motifs(os.path.join(d, s, "summary.tsv"))
            rec[f"{s}_n_motifs"] = n
            rec[f"{s}_top_consensus"] = cons
            rec[f"{s}_top_width"] = w
            rec[f"{s}_top_Evalue"] = e
            rec[f"{s}_palindromy"] = palindromic(cons)

    elif axis == "chexmix":
        for s in strains:
            rec[f"{s}_events"] = n_events(os.path.join(d, f"{s}.events"))

    else:
        for s in strains:
            c = unique_gene_counts(os.path.join(d, f"{s}.direct_regulon.tsv"))
            dp = c.get("DIRECT_PERBOX", 0)
            di = c.get("DIRECT", 0)
            rec[f"{s}_DIRECT_PERBOX"] = dp
            rec[f"{s}_DIRECT"] = di
            rec[f"{s}_direct_total"] = dp + di
            rec[f"{s}_BOUND_NOT_DE"] = c.get("BOUND_NOT_DE", 0)
            rec[f"{s}_INDIRECT"] = c.get("INDIRECT", 0)
            rec[f"{s}_genes_touched"] = c.get("_total_genes", 0)
            de = dp + di + c.get("INDIRECT", 0)
            rec[f"{s}_pct_DE_indirect"] = (
                f"{100.0*c.get('INDIRECT',0)/de:.1f}" if de else "")
        cs = crossstrain(os.path.join(d, "crossstrain_summary.tsv"))
        for k in ("BOUND_BOTH", "SHARED_DIRECT", "SHARED_BOUND",
                  "DIVERGENT_BINDING", "SHARED_INDIRECT"):
            rec[k] = cs.get(k, "")
        if axis == "orthology":
            rec["RBH_1to1"] = rbh_count(os.path.join(d, "orthologs_rbh.tsv"))

    rows.append(rec)

cols = ["axis", "combo"] + param_keys
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)

with open(out, "w") as ob:
    ob.write("\t".join(cols) + "\n")
    for r in rows:
        ob.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def spread(key):
    v = [r[key] for r in rows if str(r.get(key, "")).strip() not in ("", "None")]
    try:
        v = [float(x) for x in v]
    except ValueError:
        return None
    return (min(v), max(v), len(set(v)))


with open(logp, "w") as lg:
    lg.write(f"axis={axis} combinations={len(rows)}\n")
    lg.write("counts are UNIQUE GENES (class precedence applied), not "
             "(peak,gene) rows\n\n")
    for key in cols:
        if key in ("axis", "combo") or key in param_keys:
            continue
        sp = spread(key)
        if not sp:
            continue
        lo, hi, nuniq = sp
        verdict = ("INVARIANT" if nuniq == 1 else
                   "stable" if hi <= lo * 1.15 else "SENSITIVE")
        lg.write(f"  {key:28s} {lo:g} to {hi:g}   ({nuniq} distinct)  {verdict}\n")
