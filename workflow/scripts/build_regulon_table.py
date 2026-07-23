"""Join peak->gene assignments with RNA-seq to produce the DIRECT REGULON table.

Classification per gene:
  DIRECT           bound (peak at its promoter / TU) AND differentially expressed
  DIRECT_PERBOX    as DIRECT, and a Per-box match lies under the peak
  BOUND_NOT_DE     bound but expression unchanged (PerR present, no measured effect)
  INDIRECT         differentially expressed but NO peak assigned -> downstream effect

Also emits a per-strain summary and a strict join-rate check: if the RNA-seq
gene_ids do not match the GFF locus tags, everything downstream is silently
wrong, so the join rate is reported and a hard failure is raised at 0%.
"""
import os

assign = snakemake.input.assign         # noqa: F821
deg    = snakemake.input.deg            # noqa: F821
out    = snakemake.output.tsv           # noqa: F821
summ   = snakemake.output.summary       # noqa: F821
strain = snakemake.params.strain        # noqa: F821
logp   = snakemake.log[0]               # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return hdr, [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]

ahdr, arows = read_tsv(assign)
dhdr, drows = read_tsv(deg)
ai = {c: i for i, c in enumerate(ahdr)}
di = {c: i for i, c in enumerate(dhdr)}

rna = {}
for r in drows:
    if len(r) <= max(di.values()):
        r = r + [""] * (max(di.values()) + 1 - len(r))
    rna[r[di["gene_id"]].strip()] = r

# ---- join rate check ------------------------------------------------------
bound_genes = {r[ai["gene_id"]] for r in arows}
joined = bound_genes & set(rna)
rate = (len(joined) / len(bound_genes)) if bound_genes else 0.0

cols = ahdr + ["rna_baseMean", "rna_log2FC", "rna_padj", "deg_status",
               "in_lab_deg_list", "classification"]

counts = {}
with open(out, "w") as ob:
    ob.write("\t".join(cols) + "\n")
    for r in arows:
        g = r[ai["gene_id"]]
        rr = rna.get(g)
        def val(k):
            return rr[di[k]] if (rr and k in di and di[k] < len(rr)) else ""
        status = val("deg_status") or "no_rna"
        perbox = r[ai["perbox_in_window"]] == "True"
        if status in ("up", "down"):
            cls = "DIRECT_PERBOX" if perbox else "DIRECT"
        elif status == "ns":
            cls = "BOUND_NOT_DE"
        else:
            cls = "BOUND_NO_RNA"
        counts[cls] = counts.get(cls, 0) + 1
        ob.write("\t".join(r + [val("baseMean"), val("log2FC"), val("padj"),
                                status, val("in_lab_deg_list"), cls]) + "\n")

    # INDIRECT: DE genes with no peak assigned
    n_indirect = 0
    for g, rr in rna.items():
        if g in bound_genes:
            continue
        st = rr[di["deg_status"]] if "deg_status" in di and di["deg_status"] < len(rr) else ""
        if st in ("up", "down"):
            n_indirect += 1
            blank = [strain, "", "", "", "", "", "", "", "", g, "", "", "", "none", ""]
            ob.write("\t".join(blank + [
                rr[di["baseMean"]] if "baseMean" in di else "",
                rr[di["log2FC"]] if "log2FC" in di else "",
                rr[di["padj"]] if "padj" in di else "",
                st,
                rr[di["in_lab_deg_list"]] if "in_lab_deg_list" in di else "",
                "INDIRECT"]) + "\n")
counts["INDIRECT"] = n_indirect

with open(summ, "w") as sb:
    sb.write("strain\tcategory\tn\n")
    for k in sorted(counts):
        sb.write(f"{strain}\t{k}\t{counts[k]}\n")
    sb.write(f"{strain}\tgene_id_join_rate\t{rate:.3f}\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain} bound_genes={len(bound_genes)} "
             f"joined_to_rna={len(joined)} join_rate={rate:.3f}\n")
    for k in sorted(counts):
        lg.write(f"  {k}: {counts[k]}\n")

if bound_genes and rate == 0.0:
    raise SystemExit(
        "FATAL: no RNA-seq gene_id matched a GFF locus_tag. "
        "Check that the DEG sheet uses SAK_RS*/ID870_RS* locus tags.")