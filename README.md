# PerR ChIP-exo — reproducible analysis workflow (Group B *Streptococcus*)

A Snakemake workflow that takes the **exact zip delivered by the sequencing
core** and, in a single run, creates all directories, unpacks the archive,
auto-discovers the samples, and produces a full ChIP-exo quality-control report.
Peak calling, motif discovery, and RNA-seq integration are added in later
stages; **this release covers ingestion → QC → report**.

Biological question: does the peroxide-responsive Fur-family repressor **PerR**
directly bind the promoters of the genes it appears to regulate in GBS, and how
does that direct regulon compare between serotype Ia (**A909**) and serotype V
(**CJB111**)?

---

## What you get from one command

```
results/
├── 00_raw/                     # extracted delivery (bam, fastq.gz, bw, report)
├── 01_qc/
│   ├── header_check/           # genome-build verification per sample
│   ├── samtools/               # flagstat / idxstats / stats
│   ├── fastqc/                 # raw-read QC (R1 & R2)
│   ├── chipexoqual/            # ARC/URC, FSR, beta1/beta2 per sample
│   └── deeptools/              # fingerprint, correlation, fragment size / strain
└── 02_report/
    ├── multiqc_report.html     # aggregated FastQC + samtools
    ├── qc_scorecard.tsv        # one row per sample, with advisory flags
    └── qc_scorecard.md         # same, human-readable
```

See `docs/QC_METRICS.md` for what every metric means and its citation.

## The only thing a user edits

`config/config.yaml` → `input_zip:` (path to the core's `.zip`). Everything else
(sample identities, IP vs control, strain → genome accession) is parsed from the
CEGR file-name convention in the archive's manifest — no sample sheet to write.

Design encoded in the config, straight from the core report:
`A909 → NC_007432.1`, `CJB111 → NZ_CP063198.2`; 3 × 6×His PerR IP + 1 no-tag
control per strain (8 samples).

---

## Quick start (local / a single interactive node)

```bash
# 1. get the code
git clone https://github.com/<you>/perr-chipexo.git
cd perr-chipexo

# 2. create the driver environment (Snakemake 8 + Slurm plugin)
mamba create -n snakemake -c conda-forge -c bioconda \
    snakemake=8 python=3.12 snakemake-executor-plugin-slurm
mamba activate snakemake

# 3. point the workflow at the delivery
mkdir -p data && cp /path/to/844_LBurcham_UTenn.zip data/
#   (or edit input_zip: in config/config.yaml)

# 4. dry-run to see the DAG, then run
snakemake -n --sdm conda                      # dry run
snakemake --sdm conda --cores 8               # execute locally
```

`--sdm conda` (software-deployment-method) makes Snakemake build each rule's
pinned conda env from `workflow/envs/`, so the only thing you install by hand is
Snakemake itself.

---

## Running on ISAAC-NG (UTK, Slurm)

ISAAC-NG uses Slurm; this repo ships a workflow profile at `profiles/isaac/`.

```bash
# on a login node
module load anaconda3                          # or your miniconda module
mamba create -n snakemake -c conda-forge -c bioconda \
    snakemake=8 python=3.12 snakemake-executor-plugin-slurm
mamba activate snakemake

# edit profiles/isaac/config.yaml -> set slurm_account / slurm_partition / qos
#   (defaults mirror your history: acf-utk0009 / condo-rmccord / condo)

# build conda envs ONCE on the login node (needs internet); compute nodes reuse
snakemake --sdm conda --conda-create-envs-only --cores 1

# launch: Snakemake submits each job to Slurm for you
snakemake --workflow-profile profiles/isaac
```

Notes for ISAAC:
- Run the workflow from a persistent path (e.g. your project/scratch space),
  not `$HOME` quota-limited dirs; the extracted delivery is ~36 GB.
- The profile submits jobs on your behalf — you do **not** wrap the call in
  `sbatch`. Run it inside `tmux`/`screen` on the login node (or as a single
  lightweight batch job that itself runs `snakemake`).
- Conda envs are built on the login node (internet); compute nodes then reuse
  them offline. The QC stage needs **no** internet on compute nodes.

---

## Repository layout

```
perr-chipexo/
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── config.yaml             # <-- the one file you edit
├── workflow/
│   ├── Snakefile               # targets + includes
│   ├── rules/
│   │   ├── common.smk          # sample discovery from the zip manifest
│   │   ├── setup.smk           # extraction (zip-bomb-safe)
│   │   ├── qc.smk              # samtools / fastqc / chipexoqual / deeptools
│   │   └── report.smk          # multiqc + scorecard
│   ├── envs/                   # per-rule pinned conda envs
│   └── scripts/                # python + R helper scripts
├── profiles/
│   └── isaac/config.yaml       # Snakemake 8 Slurm profile for ISAAC-NG
└── docs/
    └── QC_METRICS.md           # metric definitions + citations
```

---

## Acknowledgements / citations to include in the paper

- Sequencing performed at the **Cornell BRC Epigenomics Facility**
  (RRID:SCR_021287).
- **ScriptManager** (Lang et al. 2022; RRID:SCR_021797) for downstream
  composite/heatmap visualisation (later stage).
- **ChIPexoQual** — Welch et al., *Nucleic Acids Research* 2017.
- **deepTools** — Ramírez et al., *Nucleic Acids Research* 2016.
- Peak calling (later stage): **ChExMix** (Yamada et al. 2019) with
  cross-checks; **DEOCSU** deep-learning caller (ChIP-mini, *NAR* 2025).

## Roadmap (next stages)
1. Reference download + annotation (`reference.download: true`).
2. Peak calling: ChExMix (primary) + MACE/GEM + DEOCSU concordance.
3. Motif: MEME-ChIP vs canonical 7-1-7 PerR box; FIMO genome scan.
4. Visualisation: IGV loci + ScriptManager/deepTools metaplots.
5. ChIP-exo × RNA-seq direct-target table; cross-serotype comparison.
