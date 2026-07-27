"""Join peak->gene assignments with RNA-seq to produce the regulon table.

Classification per gene:
  DIRECT           bound (peak at its promoter / TU) AND differentially expressed
  DIRECT_PERBOX    as DIRECT, and a Per-box match lies under the peak
  BOUND_NOT_DE     bound but expression unchanged (PerR present, no measured effect)
  BOUND_NO_RNA     bound, but no RNA-seq measurement for this gene
  INDIRECT         differentially expressed but NO peak assigned -> downstream effect

Also emits a per-strain summary and a strict join-rate check: if the RNA-seq
gene_ids do not match the GFF locus tags, everything downstream is silently
wrong, so the join rate is reported and a hard failure is raised at 0%.

-----------------------------------------------------------------------------
v2 CHANGES
-----------------------------------------------------------------------------
1. THE GFF IS NOW AN INPUT.  In v1 the INDIRECT rows were written from a
   hard-coded blank list, so every DE-but-unbound gene carried a locus tag and
   NOTHING else -- no gene name, no product, no strand.  That was 284 of the
   487 rows in the cross-strain table (58%), including the entire set backing
   the "~76% of DEGs are indirect" result.  A supplementary table of bare locus
   tags is not usable.  INDIRECT rows are now annotated from the GFF exactly
   like bound rows.

2. NEW TRAILING COLUMN `gene_biotype` (protein_coding / rRNA / tRNA / ...).
   Non-coding genes have no CDS, therefore no protein in the .faa, therefore
   they can NEVER have a reciprocal best hit no matter how conserved they are.
   Without this column they are indistinguishable downstream from genes that
   are genuinely absent from the other genome, and "NO_ORTHOLOG" conflates a
   real accessory-genome signal with an artifact of a protein-based method.
   The column is appended LAST so that any positional indexing of the existing
   columns by downstream scripts is unaffected.

3. INDIRECT rows are now built BY COLUMN NAME rather than as a fixed-length
   literal list, so adding or reordering assignment columns upstream can no
   longer silently shift fields into the wrong headers.

4. Which Per box drives DIRECT_PERBOX is configurable (`gas` | `gbs`).  Until
   the GBS motif scan is wired in, `gbs` falls back to the GAS column and says
   so in the log, so this script is safe to install before that change lands.
"""
import os
import re

assign = snakemake.input.assign          # noqa: F821
deg = snakemake.input.deg                # noqa: F821
gff = snakemake.input.gff                # noqa: F821
out = snakemake.output.tsv               # noqa: F821
summ = snakemake.output.summary          # noqa: F821
strain = snakemake.params.strain         # noqa: F821
perbox_source = str(getattr(snakemake.params, "perbox_source", "gas")).lower()  # noqa: F821
logp = snakemake.log[0]                  # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return hdr, [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]


def parse_gff(path):
    """locus_tag -> {gene_name, product, strand, biotype} for EVERY gene,
    coding or not.  `gene` features carry gene= and gene_biotype=; the child
    feature (CDS / rRNA / tRNA / ncRNA / tmRNA) carries product=."""
    info = {}

    def rec(lt):
        return info.setdefault(
            lt, {"gene_name": "", "product": "", "strand": "", "biotype": ""})

    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            lt = re.search(r"locus_tag=([^;]+)", f[8])
            if not lt:
                continue
            r = rec(lt.group(1))
            if f[2] == "gene":
                nm = re.search(r"(?:^|;)gene=([^;]+)", f[8])
                bt = re.search(r"gene_biotype=([^;]+)", f[8])
                if nm:
                    r["gene_name"] = nm.group(1)
                if bt:
                    r["biotype"] = bt.group(1)
                r["strand"] = f[6]
            elif f[2] in ("CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "RNase_P_RNA",
                          "SRP_RNA", "antisense_RNA", "riboswitch"):
                pr = re.search(r"product=([^;]+)", f[8])
                if pr and not r["product"]:
                    r["product"] = pr.group(1)
    return info


ahdr, arows = read_tsv(assign)
dhdr, drows = read_tsv(deg)
ai = {c: i for i, c in enumerate(ahdr)}
di = {c: i for i, c in enumerate(dhdr)}
ginfo = parse_gff(gff)

rna = {}
for r in drows:
    if len(r) <= max(di.values()):
        r = r + [""] * (max(di.values()) + 1 - len(r))
    rna[r[di["gene_id"]].strip()] = r

# ---- which Per box drives DIRECT_PERBOX -----------------------------------
GBS_COL, GAS_COL = "gbs_perbox_in_window", "perbox_in_window"
if perbox_source == "gbs" and GBS_COL in ai:
    perbox_col, perbox_note = GBS_COL, "GBS de novo motif"
elif perbox_source == "gbs":
    perbox_col, perbox_note = GAS_COL, "GAS box (GBS column not present yet -- fell back)"
else:
    perbox_col, perbox_note = GAS_COL, "GAS box (reference yardstick)"

# ---- join rate check ------------------------------------------------------
bound_genes = {r[ai["gene_id"]] for r in arows}
joined = bound_genes & set(rna)
rate = (len(joined) / len(bound_genes)) if bound_genes else 0.0

cols = ahdr + ["rna_baseMean", "rna_log2FC", "rna_padj", "deg_status",
               "in_lab_deg_list", "classification", "gene_biotype"]

counts, biotypes = {}, {}
n_annotated_indirect = 0

with open(out, "w") as ob:
    ob.write("\t".join(cols) + "\n")

    # ---- bound genes ------------------------------------------------------
    for r in arows:
        g = r[ai["gene_id"]]
        rr = rna.get(g)

        def val(k):
            return rr[di[k]] if (rr and k in di and di[k] < len(rr)) else ""

        status = val("deg_status") or "no_rna"
        perbox = (r[ai[perbox_col]] == "True") if perbox_col in ai else False
        if status in ("up", "down"):
            cls = "DIRECT_PERBOX" if perbox else "DIRECT"
        elif status == "ns":
            cls = "BOUND_NOT_DE"
        else:
            cls = "BOUND_NO_RNA"
        bt = ginfo.get(g, {}).get("biotype", "")
        counts[cls] = counts.get(cls, 0) + 1
        biotypes[bt or "unknown"] = biotypes.get(bt or "unknown", 0) + 1
        ob.write("\t".join(r + [val("baseMean"), val("log2FC"), val("padj"),
                                status, val("in_lab_deg_list"), cls, bt]) + "\n")

    # ---- INDIRECT: DE genes with no peak assigned -------------------------
    n_indirect = 0
    for g, rr in rna.items():
        if g in bound_genes:
            continue
        st = rr[di["deg_status"]] if ("deg_status" in di and di["deg_status"] < len(rr)) else ""
        if st not in ("up", "down"):
            continue
        n_indirect += 1
        gi = ginfo.get(g, {})
        if gi.get("gene_name") or gi.get("product"):
            n_annotated_indirect += 1

        # build by COLUMN NAME, never by fixed position
        row = [""] * len(ahdr)
        row[ai["strain"]] = strain
        row[ai["gene_id"]] = g
        row[ai["gene_name"]] = gi.get("gene_name", "")
        row[ai["product"]] = gi.get("product", "")
        row[ai["gene_strand"]] = gi.get("strand", "")
        row[ai["relationship"]] = "none"

        bt = gi.get("biotype", "")
        biotypes[bt or "unknown"] = biotypes.get(bt or "unknown", 0) + 1
        ob.write("\t".join(row + [
            rr[di["baseMean"]] if "baseMean" in di else "",
            rr[di["log2FC"]] if "log2FC" in di else "",
            rr[di["padj"]] if "padj" in di else "",
            st,
            rr[di["in_lab_deg_list"]] if "in_lab_deg_list" in di else "",
            "INDIRECT", bt]) + "\n")

counts["INDIRECT"] = n_indirect

with open(summ, "w") as sb:
    sb.write("strain\tcategory\tn\n")
    for k in sorted(counts):
        sb.write(f"{strain}\t{k}\t{counts[k]}\n")
    for k in sorted(biotypes):
        sb.write(f"{strain}\tbiotype_{k}\t{biotypes[k]}\n")
    sb.write(f"{strain}\tgene_id_join_rate\t{rate:.3f}\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain} bound_genes={len(bound_genes)} "
             f"joined_to_rna={len(joined)} join_rate={rate:.3f}\n")
    lg.write(f"DIRECT_PERBOX driven by: {perbox_note} (column '{perbox_col}')\n")
    lg.write(f"genes in GFF: {len(ginfo)}\n")
    lg.write(f"INDIRECT rows: {n_indirect}  annotated from GFF: {n_annotated_indirect}\n")
    if n_indirect and n_annotated_indirect < n_indirect:
        lg.write(f"  NOTE: {n_indirect - n_annotated_indirect} INDIRECT genes had no "
                 f"gene name or product in the GFF (hypothetical / no annotation)\n")
    for k in sorted(counts):
        lg.write(f"  {k}: {counts[k]}\n")
    for k in sorted(biotypes):
        lg.write(f"  biotype {k}: {biotypes[k]}\n")

if bound_genes and rate == 0.0:
    raise SystemExit(
        "FATAL: no RNA-seq gene_id matched a GFF locus_tag. "
        "Check that the DEG sheet uses SAK_RS*/ID870_RS* locus tags.")