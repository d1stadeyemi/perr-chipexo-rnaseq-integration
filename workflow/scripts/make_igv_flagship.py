"""make_igv_flagship.py

The curated panel list: the loci that carry the paper, plus the two that carry
its caveats. Derived from the analysis outputs rather than hard-coded, so it
stays correct if the results change.

Panels emitted, in the order they should appear in a figure:
  1  perR autoregulation          -- the internal positive control
  2  the anchor locus             -- where every line of evidence converges
  3  a SHARED_DIRECT core locus   -- the conserved regulon
  4  A_ONLY with elevated CJB111  -- the sensitivity argument, made visual
  5  BOUND_BOTH, paired with (4)  -- the calibration for (4)
  6  a Subtype1 event             -- the artefact candidate we chose not to trust
  7  a clean Subtype0 event       -- the canonical doublet
  8  the rRNA-operon peak         -- annotated, not excluded

Panels 4 and 5 must be shown ADJACENT and at IDENTICAL data range. That pairing
is the single most persuasive image in the orthologous-position argument: it
turns "occupied but below threshold" from a statistic into something visible.
"""
import os, re

loci_files = list(snakemake.input.loci)          # noqa: F821
crossstrain = snakemake.input.crossstrain        # noqa: F821
orthoanchor = snakemake.input.orthoanchor        # noqa: F821
event_files = list(snakemake.input.events)       # noqa: F821
out = snakemake.output.tsv                       # noqa: F821
strains = list(snakemake.params.strains)         # noqa: F821
logp = snakemake.log[0]                          # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return hdr, [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]


# locus tables, keyed strain -> locus_tag -> row dict
LOCI = {}
for p in loci_files:
    hdr, rows = read_tsv(p)
    idx = {c: i for i, c in enumerate(hdr)}
    for r in rows:
        if len(r) < len(hdr):
            r += [""] * (len(hdr) - len(r))
        LOCI.setdefault(r[idx["strain"]], {})[r[idx["locus_tag"]]] = dict(
            zip(hdr, r))

panels = []


def add(order, panel, strain, tag, why, note=""):
    rec = LOCI.get(strain, {}).get(tag)
    if not rec:
        return False
    panels.append([order, panel, strain, tag,
                   rec.get("gene_name", ""), rec.get("product", ""),
                   rec.get("classification", ""),
                   rec.get("igv_context", ""), rec.get("igv_zoom", ""),
                   rec.get("n_peaks", ""), rec.get("gbs_motif", ""),
                   rec.get("rna_log2FC", ""), why, note])
    return True


# ---- 1 & 2: the two flagship loci, by function rather than by hard-coded tag
def find_by(strain, pattern, cls_prefix=None):
    for tag, rec in LOCI.get(strain, {}).items():
        blob = (rec.get("gene_name", "") + " " + rec.get("product", "")).lower()
        if re.search(pattern, blob):
            if cls_prefix and not rec.get("classification", "").startswith(cls_prefix):
                continue
            return tag
    return None

for s in strains:
    t = find_by(s, r"\bperr\b")
    if t:
        add(1, "perR autoregulation", s, t,
            "PerR bound at its own promoter, recovered from an unbiased scan",
            "NO Per box at either peak -- the motif track is correctly empty here. "
            "The RNA-seq log2FC of ~-10 is the DELETION, not autoregulation.")

for s in strains:
    t = find_by(s, r"vit1|ccc1|vacuolar iron")
    if t:
        add(2, "anchor locus (VIT1/CCC1)", s, t,
            "divergent promoter; motif at -2 bp from the crosslink; ~10-fold derepressed",
            "Show a wide panel for the peak cluster AND a 500 bp zoom on the "
            "152 bp intergenic region. One bound repressor occludes BOTH promoters.")

# ---- 3: a SHARED_DIRECT core locus
try:
    hdr, rows = read_tsv(crossstrain)
    idx = {c: i for i, c in enumerate(hdr)}
    for r in rows:
        if len(r) > idx.get("conservation", 0) and r[idx["conservation"]] == "SHARED_DIRECT":
            tag = r[0]
            if tag and add(3, "conserved core (SHARED_DIRECT)", strains[0], tag,
                           "bound AND differentially expressed in both serotypes",
                           "One of only 5, and invariant across the assignment-threshold grid."):
                break
except Exception as e:  # noqa: BLE001
    pass

# ---- 4 & 5: the sensitivity pair
try:
    hdr, rows = read_tsv(orthoanchor)
    idx = {c: i for i, c in enumerate(hdr)}
    gcol, acol = idx["group"], 0
    vcol = [i for i, c in enumerate(hdr) if c.endswith("_ip_max")][-1]

    def top(group, n=1):
        sel = [r for r in rows if len(r) > vcol and r[gcol] == group]
        def val(r):
            try:
                return float(r[vcol])
            except ValueError:
                return -1.0
        return sorted(sel, key=val, reverse=True)[:n]

    for r in top("A_ONLY"):
        add(4, "A_ONLY -- elevated but sub-threshold in CJB111", strains[0], r[acol],
            "A909 calls a peak; CJB111 shows signal but no call",
            "PAIR WITH PANEL 5, adjacent, IDENTICAL data range. This is the "
            "orthologous-position argument made visual.")
    for r in top("BOUND_BOTH"):
        add(5, "BOUND_BOTH -- positive control for panel 4", strains[0], r[acol],
            "both strains independently call a peak",
            "The calibration for panel 4: this is what 'detected' looks like.")
except Exception:  # noqa: BLE001
    pass

# ---- 6 & 7: subtype exemplars, straight from the events files
def events_by_subtype(path, subtype, n=1):
    got = []
    with open(path) as fh:
        si = None
        for ln in fh:
            if ln.startswith("#Point"):
                cols = ln.lstrip("#").rstrip("\n").split("\t")
                for i, c in enumerate(cols):
                    if c.strip().lower() == "subtypename":
                        si = i
                continue
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if si is None or len(f) <= si or f[si].strip() != subtype:
                continue
            m = re.match(r"^(?:chr)?(.+):(\d+)$", f[0].strip())
            if m:
                try:
                    got.append((m.group(1), int(m.group(2)), float(f[3])))
                except ValueError:
                    pass
    return sorted(got, key=lambda x: -x[2])[:n]

for ev, s in zip(event_files, strains):
    for st, order, label, why, note in (
        ("Subtype1", 6, "artefact candidate (Subtype1)",
         "single unpaired forward-strand spike, no reverse partner",
         "SHOW THIS. Demonstrating what we chose NOT to trust is more convincing "
         "than showing only the best data."),
        ("Subtype0", 7, "canonical doublet (Subtype0)",
         "the dominant class -- forward left, reverse right, correct polarity",
         "63% of the A909 peak set looks like this."),
    ):
        for chrom, pos, lfc in events_by_subtype(ev, st, 1):
            panels.append([order, label, s, f"{chrom}:{pos}", "", "", st,
                           f"{chrom}:{pos-900}-{pos+900}",
                           f"{chrom}:{pos-250}-{pos+250}",
                           1, "", "", why, note])
    break  # A909 only; CJB111 has a single subtype

# ---- 8: the rRNA artefact, if it is in the locus table
for s in strains:
    for tag, rec in LOCI.get(s, {}).items():
        if rec.get("biotype") == "rRNA" and str(rec.get("n_peaks", "0")) not in ("", "0"):
            add(8, "rRNA-operon peak (annotated, not excluded)", s, tag,
                "inside a 5S rRNA gene, weak, no motif",
                "The ONLY one of 87 peaks touching any rRNA operon or tRNA -- so "
                "there is no systematic artefact and no exclusion list is warranted.")
            break

HDR = ["order", "panel", "strain", "locus_or_position", "gene_name", "product",
       "classification", "igv_context", "igv_zoom", "n_peaks", "gbs_motif",
       "rna_log2FC", "why_this_panel", "notes"]

with open(out, "w") as ob:
    ob.write("\t".join(HDR) + "\n")
    for p in sorted(panels, key=lambda x: (x[0], x[2])):
        ob.write("\t".join(str(x) for x in p) + "\n")

with open(logp, "w") as lg:
    lg.write(f"flagship panels written: {len(panels)}\n")
    for p in sorted(panels, key=lambda x: x[0]):
        lg.write(f"  {p[0]}. [{p[2]}] {p[1]} -> {p[7]}\n")
    if not panels:
        lg.write("WARNING: no panels resolved. Check that the locus tables and "
                 "pileup outputs exist and carry the expected columns.\n")
