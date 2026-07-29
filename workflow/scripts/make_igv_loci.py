"""make_igv_loci.py

Build a lookup table answering: "I found an interesting gene -- where do I look
in IGV, and what should I expect to see there?"

One row per gene touched by the analysis, carrying a ready-to-paste IGV
coordinate string plus everything already known about that locus: which peaks
were assigned to it, their fold-enrichment, whether a motif sits under them, and
the expression response.

Two coordinate strings are given per gene:
  igv_context  gene body +/- flank_context bp  -- shows flanking genes, the
               right default for a first look and for a figure panel
  igv_zoom     TSS +/- flank_zoom bp           -- resolves the strand doublet,
               which is invisible at context width

Why the TSS rather than the gene midpoint for the zoom: PerR is a repressor, so
the binding that matters is at the promoter, and for minus-strand genes that is
the gene END. Centring on the midpoint would put half the zoom panels in the
wrong place.
"""
import os
import re
from collections import defaultdict

regulon = snakemake.input.regulon            # noqa: F821
assign = snakemake.input.assign              # noqa: F821
gff = snakemake.input.gff                    # noqa: F821
out = snakemake.output.tsv                   # noqa: F821
strain = snakemake.params.strain             # noqa: F821
FLANK_C = int(snakemake.params.flank_context)   # noqa: F821
FLANK_Z = int(snakemake.params.flank_zoom)      # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


def strip_chr(n):
    return n[3:] if n.startswith("chr") else n


# ---- gene coordinates ------------------------------------------------------
genes, chrom_len = {}, {}
with open(gff) as fh:
    for ln in fh:
        if ln.startswith("##sequence-region"):
            p = ln.split()
            if len(p) >= 4:
                chrom_len[p[1]] = int(p[3])
            continue
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        lt = re.search(r"locus_tag=([^;]+)", f[8])
        if not lt:
            continue
        nm = re.search(r"(?:^|;)gene=([^;]+)", f[8])
        bt = re.search(r"gene_biotype=([^;]+)", f[8])
        genes[lt.group(1)] = dict(
            chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
            name=nm.group(1) if nm else "",
            biotype=bt.group(1) if bt else "",
        )

# ---- peaks per gene --------------------------------------------------------
peaks = defaultdict(list)
with open(assign) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ai = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(hdr):
            f += [""] * (len(hdr) - len(f))

        def av(k):
            i = ai.get(k)
            return f[i] if i is not None else ""

        peaks[av("gene_id")].append(dict(
            pos=av("peak_pos"), lfc=av("chip_log2FC"), subtype=av("subtype"),
            rel=av("relationship"), dist=av("dist_to_TSS"),
            gbs=av("gbs_perbox_in_window"), gas=av("perbox_in_window"),
            repcode=av("replication_code"),
        ))

# ---- classification and expression -----------------------------------------
best = {}
PREC = ["DIRECT_PERBOX", "DIRECT", "BOUND_NOT_DE", "BOUND_NO_RNA", "INDIRECT"]
with open(regulon) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ri = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(hdr):
            f += [""] * (len(hdr) - len(f))
        g = f[ri["gene_id"]]
        cls = f[ri["classification"]]
        rank = PREC.index(cls) if cls in PREC else 99
        if g not in best or rank < best[g][0]:
            best[g] = (rank, dict(
                cls=cls,
                product=f[ri["product"]] if "product" in ri else "",
                gene_name=f[ri["gene_name"]] if "gene_name" in ri else "",
                log2FC=f[ri["rna_log2FC"]] if "rna_log2FC" in ri else "",
                padj=f[ri["rna_padj"]] if "rna_padj" in ri else "",
                biotype=f[ri["gene_biotype"]] if "gene_biotype" in ri else "",
            ))

HDR = ["strain", "locus_tag", "gene_name", "product", "classification",
       "chrom", "start", "end", "strand", "biotype",
       "igv_context", "igv_zoom",
       "n_peaks", "peak_positions", "max_chip_log2FC", "subtypes",
       "replication_codes", "gbs_motif", "gas_motif",
       "rna_log2FC", "rna_padj"]

rows = []
for g, (_, info) in best.items():
    gc = genes.get(g)
    pk = peaks.get(g, [])

    if gc:
        chrom = strip_chr(gc["chrom"])
        clen = chrom_len.get(gc["chrom"], 10 ** 9)
        cs = max(1, gc["start"] - FLANK_C)
        ce = min(clen, gc["end"] + FLANK_C)
        tss = gc["start"] if gc["strand"] == "+" else gc["end"]
        zs = max(1, tss - FLANK_Z)
        ze = min(clen, tss + FLANK_Z)
        ctx = f"{chrom}:{cs:,}-{ce:,}".replace(",", "")
        zoom = f"{chrom}:{zs:,}-{ze:,}".replace(",", "")
    else:
        chrom = cs = ce = ctx = zoom = ""
        gc = dict(start="", end="", strand="", name="", biotype="")

    lfcs = []
    for p in pk:
        try:
            lfcs.append(float(p["lfc"]))
        except (TypeError, ValueError):
            pass

    rows.append([
        strain, g,
        info["gene_name"] or gc.get("name", ""),
        info["product"], info["cls"],
        chrom, gc.get("start", ""), gc.get("end", ""), gc.get("strand", ""),
        info["biotype"] or gc.get("biotype", ""),
        ctx, zoom,
        len(pk),
        ";".join(str(p["pos"]) for p in pk),
        f"{max(lfcs):.3f}" if lfcs else "",
        ";".join(sorted({p["subtype"] for p in pk if p["subtype"]})),
        ";".join(sorted({p["repcode"] for p in pk if p["repcode"]})),
        "True" if any(p["gbs"] == "True" for p in pk) else "False",
        "True" if any(p["gas"] == "True" for p in pk) else "False",
        info["log2FC"], info["padj"],
    ])

# bound genes first, then by absolute expression change: the useful ordering
def sort_key(r):
    bound = 0 if str(r[4]).startswith(("DIRECT", "BOUND")) else 1
    try:
        mag = -abs(float(r[19]))
    except (TypeError, ValueError):
        mag = 0.0
    return (bound, mag, r[1])


with open(out, "w") as ob:
    ob.write("\t".join(HDR) + "\n")
    for r in sorted(rows, key=sort_key):
        ob.write("\t".join(str(x) for x in r) + "\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain} rows={len(rows)}\n")
    lg.write(f"with coordinates={sum(1 for r in rows if r[5])}\n")
    lg.write(f"context flank={FLANK_C} bp  zoom flank={FLANK_Z} bp around the TSS\n")
    lg.write("USAGE: grep a locus tag or gene name, paste the igv_context or "
             "igv_zoom string into IGV's location box.\n")
