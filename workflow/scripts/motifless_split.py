"""motifless_split.py

Partition peak windows by GBS-motif status.

The split reads the SAME peak_classification.tsv that drives DIRECT_PERBOX, so
'has a motif' is defined once for the whole project. A second definition here
would be free to drift out of sync with the regulon tables, and the drift would
be silent.

Window names in the BED are chrom_position, matching the ChExMix event Point,
which is how the classification table is keyed.
"""
import os

classif = snakemake.input.classif           # noqa: F821
windows = snakemake.input.windows           # noqa: F821
out_neg = snakemake.output.bed_neg          # noqa: F821
out_pos = snakemake.output.bed_pos          # noqa: F821
out_tsv = snakemake.output.tsv              # noqa: F821
strain = snakemake.params.strain            # noqa: F821
logp = snakemake.log[0]                     # noqa: F821

for p in (out_neg, logp):
    os.makedirs(os.path.dirname(p), exist_ok=True)

status = {}
with open(classif) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    key = "gbs_perbox_in_window"
    if key not in ci:
        raise SystemExit(f"FATAL: '{key}' not in {classif}. Run the GBS FIMO "
                         "scan (rule fimo_scan_gbs) first.")
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) > ci[key]:
            status[f[0]] = (f[ci[key]] == "True")

n_pos = n_neg = n_unknown = 0
with open(windows) as fh, open(out_neg, "w") as ob_n, open(out_pos, "w") as ob_p:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        st = status.get(f[3])
        if st is None:
            n_unknown += 1
            continue
        (ob_p if st else ob_n).write(ln)
        if st:
            n_pos += 1
        else:
            n_neg += 1

total = n_pos + n_neg
with open(out_tsv, "w") as ob:
    ob.write("strain\tsubset\tn_peaks\tpct\n")
    ob.write(f"{strain}\tmotifpos\t{n_pos}\t"
             f"{100.0*n_pos/total:.1f}\n" if total else f"{strain}\tmotifpos\t0\t\n")
    ob.write(f"{strain}\tmotifless\t{n_neg}\t"
             f"{100.0*n_neg/total:.1f}\n" if total else f"{strain}\tmotifless\t0\t\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain}\n")
    lg.write(f"windows classified: {total}\n")
    lg.write(f"  motif-positive: {n_pos}\n")
    lg.write(f"  motif-less    : {n_neg}\n")
    if n_unknown:
        lg.write(f"WARNING: {n_unknown} windows had no classification row. "
                 "The window BED and the classification table may be out of "
                 "sync -- rerun classify_peaks.\n")