"""Compare the PerR regulon between two strains through 1:1 orthologs.

Each strain's regulon table is collapsed to one row per gene (class precedence
DIRECT_PERBOX > DIRECT > BOUND_NOT_DE > BOUND_NO_RNA > INDIRECT), then joined
on the RBH ortholog pairs.

Per ortholog pair we report both strains' classification and expression, and a
CONSERVATION call:
  SHARED_DIRECT      direct target in BOTH strains        (conserved core)
  A_ONLY_DIRECT      direct in strain A only
  B_ONLY_DIRECT      direct in strain B only
  SHARED_BOUND       bound in both, direct in neither
  DIVERGENT_BINDING  bound in one strain only
  SHARED_INDIRECT    differentially expressed in both, bound in neither
  OTHER              anything else

IMPORTANT: strain B (CJB111) has ~4x fewer peaks, largely because it was
prepared with half the bead/antibody slurry and sequenced shallower. So
"A_ONLY_DIRECT" must NOT be read as biology without the depth-matched control;
it is at least partly a sensitivity difference.
"""
import os

reg_a  = snakemake.input.reg_a          # noqa: F821
reg_b  = snakemake.input.reg_b          # noqa: F821
orth   = snakemake.input.orthologs      # noqa: F821
out    = snakemake.output.tsv           # noqa: F821
summ   = snakemake.output.summary       # noqa: F821
sa     = snakemake.params.strain_a      # noqa: F821
sb     = snakemake.params.strain_b      # noqa: F821
logp   = snakemake.log[0]               # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

PREC = ["DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "BOUND_NO_RNA", "INDIRECT"]

def gene_level(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        best = {}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                f += [""] * (len(hdr) - len(f))
            g = f[idx["gene_id"]]
            cls = f[idx["classification"]]
            r = PREC.index(cls) if cls in PREC else 99
            if g not in best or r < best[g][0]:
                best[g] = (r, f, idx)
    return {g: (v[1], v[2]) for g, v in best.items()}

A, B = gene_level(reg_a), gene_level(reg_b)

def get(rec, key):
    if rec is None:
        return ""
    f, idx = rec
    return f[idx[key]] if key in idx else ""

def is_direct(c): return c.startswith("DIRECT")
def is_bound(c):  return c.startswith("DIRECT") or c.startswith("BOUND")

pairs, counts = [], {}
with open(orth) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 7 or f[6] != "1:1":
            continue
        ga, gb = f[0], f[1]
        ra, rb = A.get(ga), B.get(gb)
        ca, cb = get(ra, "classification"), get(rb, "classification")
        if is_direct(ca) and is_direct(cb):        cons = "SHARED_DIRECT"
        elif is_direct(ca):                         cons = f"{sa}_ONLY_DIRECT"
        elif is_direct(cb):                         cons = f"{sb}_ONLY_DIRECT"
        elif is_bound(ca) and is_bound(cb):         cons = "SHARED_BOUND"
        elif is_bound(ca) or is_bound(cb):          cons = "DIVERGENT_BINDING"
        elif ca == "INDIRECT" and cb == "INDIRECT": cons = "SHARED_INDIRECT"
        elif ca or cb:                              cons = "OTHER"
        else:                                       continue
        counts[cons] = counts.get(cons, 0) + 1
        pairs.append([ga, gb,
                      get(ra, "gene_name") or get(rb, "gene_name"),
                      get(ra, "product") or get(rb, "product"),
                      ca, cb,
                      get(ra, "rna_log2FC"), get(rb, "rna_log2FC"),
                      get(ra, "perbox_in_window"), get(rb, "perbox_in_window"),
                      f[4], cons])

with open(out, "w") as ob:
    ob.write("\t".join([f"{sa}_gene", f"{sb}_gene", "gene_name", "product",
                        f"{sa}_class", f"{sb}_class",
                        f"{sa}_log2FC", f"{sb}_log2FC",
                        f"{sa}_perbox", f"{sb}_perbox",
                        "pident", "conservation"]) + "\n")
    for p in sorted(pairs, key=lambda x: (x[11], x[0])):
        ob.write("\t".join(str(v) for v in p) + "\n")

with open(summ, "w") as sb_:
    sb_.write("category\tn_ortholog_pairs\n")
    for k in sorted(counts):
        sb_.write(f"{k}\t{counts[k]}\n")

with open(logp, "w") as lg:
    lg.write(f"{sa} genes in regulon table: {len(A)}\n")
    lg.write(f"{sb} genes in regulon table: {len(B)}\n")
    lg.write(f"1:1 ortholog pairs compared: {len(pairs)}\n")
    for k in sorted(counts):
        lg.write(f"  {k}: {counts[k]}\n")