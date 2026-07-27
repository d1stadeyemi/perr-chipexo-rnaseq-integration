"""Compare the PerR regulon between two strains through 1:1 orthologs.

Each strain's regulon table is collapsed to one row per gene (class precedence
DIRECT_PERBOX > DIRECT > BOUND_NOT_DE > BOUND_NO_RNA > INDIRECT), then joined
on the RBH ortholog pairs.

TWO INDEPENDENT CONSERVATION AXES ARE REPORTED (this is the v2 change):

  binding       BINDING ONLY -- never consults expression.
                  BOUND_BOTH / BOUND_<A> / BOUND_<B> / BOUND_NEITHER
                This is the axis that speaks to conservation of the PerR
                cistrome, and the axis the half-slurry caveat acts on.

  conservation  BINDING x EXPRESSION combined (the original column).
                  SHARED_DIRECT      direct target in BOTH strains
                  <A>_ONLY_DIRECT    direct in A only
                  <B>_ONLY_DIRECT    direct in B only
                  SHARED_BOUND       bound in both, direct in neither
                  DIVERGENT_BINDING  bound in one strain only
                  SHARED_INDIRECT    DE in both, bound in neither
                  DE_<A>_ONLY        DE in A only, bound in neither
                  DE_<B>_ONLY        DE in B only, bound in neither
                  NO_ORTHOLOG_<A>    bound/DE in A, no RBH partner at all
                  NO_ORTHOLOG_<B>    bound/DE in B, no RBH partner at all

WHY BOTH AXES: v1 reported only `conservation`, which conflates them. A locus
bound with an intact Per box in BOTH strains was labelled "<A>_ONLY_DIRECT"
purely because the other strain's log2FC fell under the DE cutoff -- e.g. the
TetR/AcrR locus (SAK_RS01645) and the queA / VIT1-CCC1 divergent promoter
(SAK_RS04610 / SAK_RS04615), which share ONE binding site between two genes.
Binding conservation was therefore systematically understated. Read `binding`
for cistrome conservation, `conservation` for regulatory-response conservation.

WHY NO_ORTHOLOG: in v1, genes with no reciprocal best hit fell out of the join
silently and were invisible in every count. Those are exactly the accessory-
genome genes -- a PerR site on serotype-specific content is a real result the
old table could not display. They are now emitted with an empty partner column.

CAVEAT carried into the output: strain B (CJB111) has ~4x fewer peaks and was
prepared with half the standard bead/antibody slurry. BOUND_<A> and
<A>_ONLY_DIRECT are therefore confounded with SENSITIVITY, not established as
biology. Read alongside results/04_downsample/ and the orthologous-position
pileup test.

COLUMN NOTE: `<strain>_gasbox` is the GAS (S. pyogenes) Per-box consensus,
which this project showed is NOT the GBS motif (TOMTOM p=0.34). It is a
reference yardstick only. `<strain>_gbsbox` carries the de novo GBS motif
(TTMTAAATAAG) and is populated once the GBS FIMO scan is added upstream; until
then it is written empty rather than omitted, so downstream readers and figure
scripts see a stable column set.
"""
import os

reg_a = snakemake.input.reg_a           # noqa: F821
reg_b = snakemake.input.reg_b           # noqa: F821
orth = snakemake.input.orthologs        # noqa: F821
out = snakemake.output.tsv              # noqa: F821
summ = snakemake.output.summary         # noqa: F821
sa = snakemake.params.strain_a          # noqa: F821
sb = snakemake.params.strain_b          # noqa: F821
logp = snakemake.log[0]                 # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

PREC = ["DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "BOUND_NO_RNA", "INDIRECT"]


def gene_level(path):
    """Collapse a per-(peak,gene) regulon table to one row per gene, keeping the
    highest-precedence classification (a gene DIRECT at one peak is never also
    counted as BOUND_NOT_DE at another)."""
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
    """Field lookup that tolerates a column being absent (e.g. the GBS motif
    column before the GBS FIMO scan is wired in)."""
    if rec is None:
        return ""
    f, idx = rec
    return f[idx[key]] if key in idx else ""


def is_direct(c):
    return c.startswith("DIRECT")


def is_bound(c):
    return c.startswith("DIRECT") or c.startswith("BOUND")


def is_de(c):
    """DE = the classifier saw a significant expression change, whether or not
    the gene is bound. BOUND_NOT_DE and BOUND_NO_RNA are deliberately excluded."""
    return c in ("DIRECT", "DIRECT_PERBOX", "INDIRECT")


def binding_call(ca, cb):
    """Binding conservation, computed WITHOUT reference to expression."""
    ba, bb = is_bound(ca), is_bound(cb)
    if ba and bb:
        return "BOUND_BOTH"
    if ba:
        return f"BOUND_{sa}"
    if bb:
        return f"BOUND_{sb}"
    return "BOUND_NEITHER"


def row(ga, gb, ra, rb, pident, binding, cons):
    return [ga, gb,
            get(ra, "gene_name") or get(rb, "gene_name"),
            get(ra, "product") or get(rb, "product"),
            get(ra, "classification"), get(rb, "classification"),
            get(ra, "rna_log2FC"), get(rb, "rna_log2FC"),
            get(ra, "perbox_in_window"), get(rb, "perbox_in_window"),
            get(ra, "gbs_perbox_in_window"), get(rb, "gbs_perbox_in_window"),
            pident, binding, cons]


# ---------------------------------------------------------------------------
# 1. paired rows (1:1 orthologs)
# ---------------------------------------------------------------------------
pairs = []
cons_counts, bind_counts = {}, {}
seen_a, seen_b = set(), set()          # genes that appear in ANY ortholog row
n_ambiguous_skipped = 0

with open(orth) as fh:
    fh.readline()                       # header
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        ga, gb = f[0], f[1]
        # every gene with a partner of ANY quality counts as "has an ortholog",
        # so ambiguous pairs are not later mislabelled NO_ORTHOLOG
        seen_a.add(ga)
        seen_b.add(gb)
        if f[6] != "1:1":
            n_ambiguous_skipped += 1
            continue
        ra, rb = A.get(ga), B.get(gb)
        ca, cb = get(ra, "classification"), get(rb, "classification")

        bind = binding_call(ca, cb)

        if is_direct(ca) and is_direct(cb):
            cons = "SHARED_DIRECT"
        elif is_direct(ca):
            cons = f"{sa}_ONLY_DIRECT"
        elif is_direct(cb):
            cons = f"{sb}_ONLY_DIRECT"
        elif is_bound(ca) and is_bound(cb):
            cons = "SHARED_BOUND"
        elif is_bound(ca) or is_bound(cb):
            cons = "DIVERGENT_BINDING"
        elif is_de(ca) and is_de(cb):
            cons = "SHARED_INDIRECT"
        elif is_de(ca):
            cons = f"DE_{sa}_ONLY"       # was the uninformative "OTHER"
        elif is_de(cb):
            cons = f"DE_{sb}_ONLY"       # was the uninformative "OTHER"
        elif ca or cb:
            cons = "OTHER"               # should now be empty; kept as a tell-tale
        else:
            continue                     # in neither strain's table -> not reported

        cons_counts[cons] = cons_counts.get(cons, 0) + 1
        bind_counts[bind] = bind_counts.get(bind, 0) + 1
        pairs.append(row(ga, gb, ra, rb, f[4], bind, cons))

# ---------------------------------------------------------------------------
# 2. genes with NO ortholog at all (invisible in v1)
# ---------------------------------------------------------------------------
n_noorth_a = n_noorth_b = 0
for g, ra in sorted(A.items()):
    if g in seen_a:
        continue
    ca = get(ra, "classification")
    if not (is_bound(ca) or is_de(ca)):
        continue
    n_noorth_a += 1
    cons = f"NO_ORTHOLOG_{sa}"
    bind = binding_call(ca, "")
    cons_counts[cons] = cons_counts.get(cons, 0) + 1
    pairs.append(row(g, "", ra, None, "", bind, cons))

for g, rb in sorted(B.items()):
    if g in seen_b:
        continue
    cb = get(rb, "classification")
    if not (is_bound(cb) or is_de(cb)):
        continue
    n_noorth_b += 1
    cons = f"NO_ORTHOLOG_{sb}"
    bind = binding_call("", cb)
    cons_counts[cons] = cons_counts.get(cons, 0) + 1
    pairs.append(row("", g, None, rb, "", bind, cons))

# ---------------------------------------------------------------------------
# 3. write
# ---------------------------------------------------------------------------
HDR = [f"{sa}_gene", f"{sb}_gene", "gene_name", "product",
       f"{sa}_class", f"{sb}_class",
       f"{sa}_log2FC", f"{sb}_log2FC",
       f"{sa}_gasbox", f"{sb}_gasbox",
       f"{sa}_gbsbox", f"{sb}_gbsbox",
       "pident", "binding", "conservation"]

with open(out, "w") as ob:
    ob.write("\t".join(HDR) + "\n")
    for p in sorted(pairs, key=lambda x: (x[14], x[0], x[1])):
        ob.write("\t".join(str(v) for v in p) + "\n")

with open(summ, "w") as sf:
    sf.write("axis\tcategory\tn\n")
    for k in sorted(bind_counts):
        sf.write(f"binding\t{k}\t{bind_counts[k]}\n")
    for k in sorted(cons_counts):
        sf.write(f"conservation\t{k}\t{cons_counts[k]}\n")

with open(logp, "w") as lg:
    lg.write(f"{sa} genes in regulon table: {len(A)}\n")
    lg.write(f"{sb} genes in regulon table: {len(B)}\n")
    lg.write(f"1:1 ortholog pairs reported: {len(pairs) - n_noorth_a - n_noorth_b}\n")
    lg.write(f"ambiguous ortholog rows skipped: {n_ambiguous_skipped}\n")
    lg.write(f"no-ortholog rows: {sa}={n_noorth_a} {sb}={n_noorth_b}\n")
    lg.write("-- binding axis (expression-independent) --\n")
    for k in sorted(bind_counts):
        lg.write(f"  {k}: {bind_counts[k]}\n")
    lg.write("-- conservation axis (binding x expression) --\n")
    for k in sorted(cons_counts):
        lg.write(f"  {k}: {cons_counts[k]}\n")