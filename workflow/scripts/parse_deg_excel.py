"""Parse the Burcham-lab DEG workbook into tidy per-strain TSVs.

The workbook has one 'raw output' sheet and one 'DEGs' sheet per strain.
We use the RAW OUTPUT sheet as the analysis input because it contains EVERY
gene with its statistics — this lets us report expression for any bound gene,
including ones that are bound but NOT differentially expressed (a biologically
meaningful category the pre-filtered DEG sheet would hide).

The 'DEGs' sheet is read only to VERIFY that our thresholds reproduce the
lab's own significant-gene list; a mismatch is reported in the log.

Column set (baseMean / log2_fold_change / lfcSE / stat / p_value / q_value) is
characteristic of DESeq2.
"""
import os
import pandas as pd

xlsx      = snakemake.input.xlsx                 # noqa: F821
raw_sheet = snakemake.params.raw_sheet           # noqa: F821
deg_sheet = snakemake.params.deg_sheet           # noqa: F821
strain    = snakemake.params.strain              # noqa: F821
padj_thr  = float(snakemake.params.padj)         # noqa: F821
lfc_thr   = float(snakemake.params.log2fc)       # noqa: F821
out_tsv   = snakemake.output.tsv                 # noqa: F821
logp      = snakemake.log[0]                     # noqa: F821

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

def norm(df):
    """Normalise column names to a canonical set."""
    ren = {}
    for c in df.columns:
        k = str(c).strip().lower().replace(" ", "_")
        if k in ("gene_id", "geneid", "gene"):            ren[c] = "gene_id"
        elif k in ("basemean", "base_mean"):              ren[c] = "baseMean"
        elif k.startswith("log2"):                        ren[c] = "log2FC"
        elif k in ("lfcse",):                             ren[c] = "lfcSE"
        elif k in ("stat",):                              ren[c] = "stat"
        elif k in ("p_value", "pvalue", "pval"):          ren[c] = "pvalue"
        elif k in ("q_value", "qvalue", "padj", "fdr"):   ren[c] = "padj"
    return df.rename(columns=ren)

with open(logp, "w") as lg:
    raw = norm(pd.read_excel(xlsx, sheet_name=raw_sheet))
    lg.write(f"strain={strain} raw_sheet='{raw_sheet}' rows={len(raw)}\n")
    lg.write(f"columns={list(raw.columns)}\n")

    need = {"gene_id", "log2FC", "padj"}
    missing = need - set(raw.columns)
    if missing:
        raise SystemExit(f"missing expected columns in '{raw_sheet}': {missing}")

    raw["gene_id"] = raw["gene_id"].astype(str).str.strip()
    for c in ("baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"):
        if c in raw.columns:
            raw[c] = pd.to_numeric(raw[c], errors="coerce")

    # our own significance call (thresholds in config)
    sig = raw["padj"].notna() & (raw["padj"] < padj_thr) & (raw["log2FC"].abs() >= lfc_thr)
    raw["deg_status"] = "ns"
    raw.loc[sig & (raw["log2FC"] > 0), "deg_status"] = "up"
    raw.loc[sig & (raw["log2FC"] < 0), "deg_status"] = "down"
    lg.write(f"our_thresholds padj<{padj_thr} |log2FC|>={lfc_thr} -> "
             f"up={(raw.deg_status=='up').sum()} down={(raw.deg_status=='down').sum()}\n")

    # cross-check against the lab's own DEG sheet
    try:
        degs = norm(pd.read_excel(xlsx, sheet_name=deg_sheet))
        lab_ids = set(degs["gene_id"].astype(str).str.strip())
        ours    = set(raw.loc[sig, "gene_id"])
        lg.write(f"lab_DEG_sheet='{deg_sheet}' n={len(lab_ids)} ours={len(ours)} "
                 f"shared={len(lab_ids & ours)} lab_only={len(lab_ids - ours)} "
                 f"ours_only={len(ours - lab_ids)}\n")
        raw["in_lab_deg_list"] = raw["gene_id"].isin(lab_ids)
    except Exception as e:                                    # noqa: BLE001
        lg.write(f"lab DEG sheet unavailable ({e}); continuing with our calls\n")
        raw["in_lab_deg_list"] = pd.NA

    raw.insert(0, "strain", strain)
    raw.to_csv(out_tsv, sep="\t", index=False)
    lg.write(f"wrote {out_tsv}\n")