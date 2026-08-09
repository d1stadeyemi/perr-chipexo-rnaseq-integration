"""Assign ChExMix binding events to genes / transcription units.

WHY NOT "NEAREST GENE":
  The strongest PerR site in this dataset (A909 874,861) sits in a 152 bp
  DIVERGENT intergenic region between queA (- strand) and SAK_RS04615
  (+ strand). A nearest-gene rule would arbitrarily pick one; a repressor
  bound there occludes BOTH promoters. So we assign a peak to EVERY gene
  whose transcription start site (TSS) lies within `max_upstream` bp
  downstream of the peak, on either strand. Divergent pairs get both.

RELATIONSHIPS reported per (peak, gene):
  upstream : peak lies 5' of the gene, within max_upstream of its TSS
             (the promoter-occlusion case; what a repressor should do)
  inside   : peak falls within the gene body
  operon   : gene is co-transcribed downstream of an 'upstream'-assigned gene
             (same strand, intergenic gap <= operon_max_gap, chained)

TSS = start for + strand genes, end for - strand genes.
"""
import os
import re

gff        = snakemake.input.gff                      # noqa: F821
events     = snakemake.input.events                   # noqa: F821
repcodes   = snakemake.input.repcodes                 # noqa: F821
# Reproducibility filter for GENE ASSIGNMENT only. Peak calling and motif
# discovery are unaffected -- discovery keeps every sequence, where statistical
# power matters, while binding claims rest only on reproducible peaks.
MIN_REPS = int(getattr(snakemake.params, "min_reps", 1))   # noqa: F821
classif    = snakemake.input.classification           # noqa: F821
out_tsv    = snakemake.output.tsv                     # noqa: F821
strain     = snakemake.params.strain                  # noqa: F821
max_up     = int(snakemake.params.max_upstream)       # noqa: F821
op_gap     = int(snakemake.params.operon_max_gap)     # noqa: F821
logp       = snakemake.log[0]                         # noqa: F821

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

# ---- genes from GFF -------------------------------------------------------
genes = []          # (chrom, start, end, strand, locus_tag, name)
products = {}       # locus_tag -> product (from CDS)
with open(gff) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom, feat, start, end, strand, attr = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        lt = re.search(r"locus_tag=([^;]+)", attr)
        if not lt:
            continue
        lt = lt.group(1)
        if feat == "gene":
            nm = re.search(r"(?:^|;)gene=([^;]+)", attr)
            genes.append((chrom, start, end, strand, lt, nm.group(1) if nm else ""))
        elif feat == "CDS":
            pr = re.search(r"product=([^;]+)", attr)
            if pr and lt not in products:
                products[lt] = pr.group(1)

genes.sort(key=lambda g: (g[0], g[1]))

def tss(g):
    return g[1] if g[3] == "+" else g[2]

# ---- operon chains: same strand, small intergenic gap ---------------------
next_in_operon = {}
for i in range(len(genes) - 1):
    a, b = genes[i], genes[i + 1]
    if a[0] == b[0] and a[3] == b[3] == "+" and (b[1] - a[2]) <= op_gap:
        next_in_operon[a[4]] = b[4]
for i in range(len(genes) - 1, 0, -1):
    a, b = genes[i], genes[i - 1]      # for - strand, transcription runs right->left
    if a[0] == b[0] and a[3] == b[3] == "-" and (a[1] - b[2]) <= op_gap:
        next_in_operon[a[4]] = b[4]

# ---- peaks ---------------------------------------------------------------
code = {}
with open(repcodes) as fh:
    for ln in fh:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip() != "X":
            code[f[0].strip()] = f[1].strip()

# Read by COLUMN NAME. v1 indexed f[5]/f[6] positionally, which silently
# breaks the moment a column is inserted upstream.
perbox = {}
with open(classif) as fh:
    chdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(chdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")

        def cv(key):
            i = ci.get(key)
            return f[i] if (i is not None and i < len(f)) else ""

        perbox[f[0]] = (cv("perbox_in_window"), cv("n_perbox_hits"),
                        cv("gbs_perbox_in_window"), cv("n_gbs_perbox_hits"))

peaks = []
with open(events) as fh:
    sub_idx = None
    for ln in fh:
        if ln.startswith("#Point"):
            cols = ln.lstrip("#").rstrip("\n").split("\t")
            for i, c in enumerate(cols):
                if c.strip().lower() == "subtypename":
                    sub_idx = i
            continue
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        point = f[0].strip()
        m = re.match(r"^(?:chr)?(.+):(\d+)$", point)
        if not m:
            continue
        chrom, pos = m.group(1), int(m.group(2))
        peaks.append({
            "point": point, "chrom": chrom, "pos": pos,
            "log2FC": f[3] if len(f) > 3 else "",
            "subtype": f[sub_idx] if (sub_idx is not None and len(f) > sub_idx) else "NA",
            "repcode": code.get(point, "NA"),
        })

# ---- assign --------------------------------------------------------------
rows = []
n_assigned = 0
n_filtered = 0
for pk in peaks:
    # Peaks called in fewer than MIN_REPS replicates are retained in the peak
    # files and used for motif discovery, but are not assigned to a gene.
    # Measured on A909, median fold enrichment rises with replication code
    # (0.80 / 0.91 / 1.30 for 1, 2 and 3 replicates), so this removes the
    # weakest tier. Unparseable codes ("NA", from peaks with no matching row in
    # the replication-codes table) are KEPT rather than silently dropped.
    try:
        if int(pk["repcode"]) < MIN_REPS:
            n_filtered += 1
            continue
    except (ValueError, TypeError, KeyError):
        pass
    hits = []
    for g in genes:
        if g[0] != pk["chrom"]:
            continue
        t = tss(g)
        if g[1] <= pk["pos"] <= g[2]:
            hits.append((g, "inside", pk["pos"] - t if g[3] == "+" else t - pk["pos"]))
            continue
        # upstream = peak on the 5' side, within max_up of the TSS
        d = (t - pk["pos"]) if g[3] == "+" else (pk["pos"] - t)
        if 0 < d <= max_up:
            hits.append((g, "upstream", -d))
    # operon extension from upstream-assigned genes
    seen = {h[0][4] for h in hits}
    for g, rel, dist in list(hits):
        if rel != "upstream":
            continue
        cur = g[4]
        while cur in next_in_operon:
            nxt = next_in_operon[cur]
            if nxt in seen:
                break
            gg = next(x for x in genes if x[4] == nxt)
            hits.append((gg, "operon", ""))
            seen.add(nxt)
            cur = nxt
    if hits:
        n_assigned += 1
    pb = perbox.get(f"{pk['chrom']}_{pk['pos']}", ("", "", "", ""))
    for g, rel, dist in hits:
        rows.append([
            strain, pk["point"], pk["chrom"], pk["pos"], pk["subtype"],
            pk["log2FC"], pk["repcode"], pb[0], pb[1],
            g[4], g[5], products.get(g[4], ""), g[3], rel, dist,
            pb[2], pb[3],
        ])

# The two GBS columns are APPENDED, not inserted, so existing positional
# indexing of this table (gene_id is still $10, relationship still $14)
# continues to work.
hdr = ["strain", "peak", "chrom", "peak_pos", "subtype", "chip_log2FC",
       "replication_code", "perbox_in_window", "n_perbox_hits",
       "gene_id", "gene_name", "product", "gene_strand",
       "relationship", "dist_to_TSS",
       "gbs_perbox_in_window", "n_gbs_perbox_hits"]
with open(out_tsv, "w") as ob:
    ob.write("\t".join(hdr) + "\n")
    for r in rows:
        ob.write("\t".join(str(x) for x in r) + "\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain} genes={len(genes)} peaks={len(peaks)} "
             f"peaks_with_gene={n_assigned} assignments={len(rows)} "
             f"max_upstream={max_up} operon_gap={op_gap}\n")
    lg.write(f"reproducibility filter: min_reps={MIN_REPS}, "
             f"{n_filtered} peaks excluded from gene assignment\n")