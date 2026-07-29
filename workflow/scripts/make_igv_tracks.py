"""Build IGV-ready annotation tracks from the analysis outputs.

THE PROBLEM THIS SOLVES FIRST
  ChExMix writes event coordinates as 'chrNC_007432.1:874861'. The reference
  FASTA, the GFF and every bigWig use 'NC_007432.1'. A track whose sequence
  names do not match the loaded genome loads SILENTLY and displays NOTHING --
  IGV does not warn. This has cost more debugging time than any other single
  issue in the project, so every track written here is prefix-normalised.

WHAT ELSE THIS ADDS
  The delivered .bed carries only coordinates. These versions put the event's
  statistics into the BED name field, so hovering a peak in IGV shows its
  fold-enrichment, ChExMix subtype and motif status without leaving the
  browser. That turns the peak track from a set of tick marks into something
  you can interrogate.

COORDINATE CONVENTIONS
  ChExMix Points are 1-based. BED is 0-based half-open. A single-base feature
  at 1-based position p is therefore written as (p-1, p). Getting this wrong
  shifts every peak by one base relative to the coverage, which is visible at
  ChIP-exo resolution.
"""
import os
import re

events = snakemake.input.events              # noqa: F821
events_ds = snakemake.input.events_ds        # noqa: F821
classif = snakemake.input.classif            # noqa: F821
fimo_gbs = snakemake.input.fimo_gbs          # noqa: F821
fimo_gas = snakemake.input.fimo_gas          # noqa: F821
regulon = snakemake.input.regulon            # noqa: F821
gff = snakemake.input.gff                    # noqa: F821

out_peaks = snakemake.output.peaks           # noqa: F821
out_peaks_ds = snakemake.output.peaks_ds     # noqa: F821
out_gbs = snakemake.output.motif_gbs         # noqa: F821
out_gas = snakemake.output.motif_gas         # noqa: F821
out_reg = snakemake.output.regulon           # noqa: F821

strain = snakemake.params.strain             # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

for p in (out_peaks, logp):
    os.makedirs(os.path.dirname(p), exist_ok=True)


def strip_chr(name):
    """Normalise a sequence name to the form used by the FASTA and GFF."""
    return name[3:] if name.startswith("chr") else name


# ---------------------------------------------------------------------------
# peak classification, keyed by the window name (chrom_position)
# ---------------------------------------------------------------------------
def read_classification(path):
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")

            def cv(k):
                i = ci.get(k)
                return f[i] if (i is not None and i < len(f)) else ""

            key = f[0]
            d[key] = {
                "subtype": cv("subtype"),
                "gbs": cv("gbs_perbox_in_window"),
                "gas": cv("perbox_in_window"),
                "offset": cv("gbs_best_offset"),
            }
    return d


CLS = read_classification(classif)


def parse_events(path):
    """-> [(chrom_no_prefix, pos_1based, log2FC, subtype)]"""
    rows = []
    with open(path) as fh:
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
            m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0].strip())
            if not m:
                continue
            chrom, pos = m.group(1), int(m.group(2))
            try:
                lfc = float(f[3])
            except (IndexError, ValueError):
                lfc = float("nan")
            st = f[sub_idx].strip() if (sub_idx is not None and len(f) > sub_idx) else "NA"
            rows.append((chrom, pos, lfc, st))
    return rows


def write_peak_bed(rows, path, tag):
    """BED6. The name field carries the statistics so IGV tooltips are useful.
    Score is scaled to BED's 0-1000 range so IGV can shade by enrichment."""
    with open(path, "w") as ob:
        ob.write(f'track name="{strain} {tag}" description="ChExMix events" '
                 f'useScore=1 itemRgb="Off"\n')
        for chrom, pos, lfc, st in rows:
            info = CLS.get(f"{chrom}_{pos}", {})
            motif = ""
            if info.get("gbs") == "True":
                off = info.get("offset", "")
                motif = f"|GBSmotif{'@' + off + 'bp' if off else ''}"
            elif info.get("gas") == "True":
                motif = "|GASbox"
            name = f"{st}|log2FC={lfc:.3f}{motif}"
            score = 1000 if lfc > 4 else max(0, min(1000, int(lfc * 240)))
            ob.write(f"{chrom}\t{pos - 1}\t{pos}\t{name}\t{score}\t.\n")
    return len(rows)


n_all = write_peak_bed(parse_events(events), out_peaks, "peaks (full depth)")
n_ds = write_peak_bed(parse_events(events_ds), out_peaks_ds, "peaks (depth-matched)")


# ---------------------------------------------------------------------------
# FIMO GFF -> prefix-normalised
# ---------------------------------------------------------------------------
def normalise_gff(src, dst, label):
    n = 0
    with open(src) as fh, open(dst, "w") as ob:
        ob.write(f'##gff-version 3\n# {strain} {label}\n')
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            f[0] = strip_chr(f[0])
            ob.write("\t".join(f) + "\n")
            n += 1
    return n


n_gbs = normalise_gff(fimo_gbs, out_gbs, "de novo GBS PerR motif (FIMO p<1e-4)")
n_gas = normalise_gff(fimo_gas, out_gas, "GAS Per box reference (FIMO p<1e-4)")


# ---------------------------------------------------------------------------
# direct-regulon BED: gene bodies of the DIRECT / DIRECT_PERBOX genes
# ---------------------------------------------------------------------------
def gene_coords(path):
    g = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            lt = re.search(r"locus_tag=([^;]+)", f[8])
            nm = re.search(r"(?:^|;)gene=([^;]+)", f[8])
            if lt:
                g[lt.group(1)] = (f[0], int(f[3]), int(f[4]), f[6],
                                  nm.group(1) if nm else "")
    return g


GC = gene_coords(gff)

direct = {}
with open(regulon) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ri = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(hdr):
            f += [""] * (len(hdr) - len(f))
        cls = f[ri["classification"]]
        if not cls.startswith("DIRECT"):
            continue
        gid = f[ri["gene_id"]]
        lfc = f[ri["rna_log2FC"]] if "rna_log2FC" in ri else ""
        prev = direct.get(gid)
        # DIRECT_PERBOX outranks DIRECT if a gene appears under both
        if prev is None or (cls == "DIRECT_PERBOX" and prev[0] != "DIRECT_PERBOX"):
            direct[gid] = (cls, lfc, f[ri["product"]] if "product" in ri else "")

with open(out_reg, "w") as ob:
    ob.write(f'track name="{strain} direct regulon" '
             f'description="DIRECT and DIRECT_PERBOX genes" itemRgb="On"\n')
    n_reg = 0
    for gid, (cls, lfc, prod) in sorted(direct.items()):
        if gid not in GC:
            continue
        chrom, s, e, strand_, nm = GC[gid]
        label = f"{gid}"
        if nm:
            label += f" ({nm})"
        label += f"|{cls}"
        if lfc:
            try:
                label += f"|log2FC={float(lfc):.2f}"
            except ValueError:
                pass
        # rust for motif-supported, teal for the rest
        rgb = "168,71,28" if cls == "DIRECT_PERBOX" else "14,110,107"
        ob.write(f"{strip_chr(chrom)}\t{s - 1}\t{e}\t{label}\t0\t{strand_}\t"
                 f"{s - 1}\t{e}\t{rgb}\n")
        n_reg += 1

with open(logp, "w") as lg:
    lg.write(f"strain={strain}\n")
    lg.write(f"peaks (full depth)    : {n_all}\n")
    lg.write(f"peaks (depth-matched) : {n_ds}\n")
    lg.write(f"GBS motif hits        : {n_gbs}\n")
    lg.write(f"GAS Per-box hits      : {n_gas}\n")
    lg.write(f"direct regulon genes  : {n_reg}\n")
    lg.write("all sequence names normalised to the FASTA/GFF convention "
             "(leading 'chr' stripped where present)\n")
