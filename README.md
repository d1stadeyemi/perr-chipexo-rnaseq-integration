# PerR regulation impacts GBS infection

This repo is a Snakemake workflow for ChIP-exo analysis of the PerR transcriptional regulator in *Streptococcus agalactiae* (Group B *Streptococcus*).

This workflow takes the raw archive delivered by the sequencing core and reproduces every result in the associated study: quality control, peak calling, de novo motif discovery, RNA-seq integration, cross-serotype orthology, and the robustness analyses that support each reported number.

This document describes the computational workflow. For the biology and the interpretation of results, see the publication linked at the bottom of this page.

---

## Overview

Two serotypes of *S. agalactiae* were profiled by ChIP-exo using a 6xHis-tagged PerR construct, with three immunoprecipitation replicates and one no-tag control each. This workflow integrates those data with matched RNA-seq from a ΔperR background to separate genes that PerR binds directly from genes that respond to its loss indirectly.

| | |
|---|---|
| Organism | *Streptococcus agalactiae*, serotypes Ia (A909) and V (CJB111) |
| Assemblies | `GCF_000012705.1`, `GCF_015221735.2` |
| Assay | ChIP-exo, paired-end, 8 samples |
| Workflow engine | Snakemake 9, per-rule conda environments |
| Input | one `.zip` from the sequencing core, plus an RNA-seq workbook |
| Runtime | approximately 6 to 9 hours on 8 cores |

Everything downstream of the archive is derived automatically. Sample identities, IP versus control assignment, replicate numbering, and the strain to genome mapping are all parsed from the file naming convention in the delivery manifest. 

---

## Repository structure

```
perr_chipexo/
├── config/
│   └── config.yaml              all parameters, fully commented
├── workflow/
│   ├── Snakefile                module includes and the `all` target
│   ├── rules/                   16 rule modules, one per analysis stage
│   ├── scripts/                 ~25 Python and R scripts called by rules
│   └── envs/                    conda environment specifications
│       └── locks/               fully resolved package lists for exact rebuilds
├── profiles/
│   └── isaac/                   SLURM profile, adapt for your cluster
├── resources/                   reference genomes and motif databases (downloaded)
├── data/                        input archive and RNA-seq workbook (user supplied)
├── results/                     all outputs, numbered by stage
├── docs/                        supplementary documentation
├── CITATIONS.md                 literature references for every tool
└── README.md                    this file
```

The `results/` tree is created by the workflow and mirrors the pipeline stages:

```
results/
├── 00_raw/            extracted delivery
├── 01_qc/             per-sample quality control
├── 02_report/         MultiQC report and QC scorecard
├── 03_peaks/          peak calls, segregated by caller
├── 04_downsample/     depth-matched re-analysis
├── 05_motif/          motif discovery, scanning, database comparison
├── 06_integration/    peak to gene assignment, regulon tables
├── 07_figures/        publication figures
├── 08_orthology/      reciprocal best hits, cross-serotype comparison
├── 09_pileup/         orthologous-position signal test
├── 10_igv/            genome browser tracks and preconfigured sessions
├── 11_sensitivity/    parameter robustness analysis
├── 12_concordance/    cross-caller agreement
├── 13_motifless/      analysis of peaks lacking the discovered motif
└── logs/              one log per rule execution
```

---

## How to reproduce

These instructions assume no prior familiarity with the project. Following them in order reproduces every result from the raw archive.

### 1. Compute requirements

| Resource | Minimum | Notes |
|---|---|---|
| CPU cores | 8 | the workflow parallelises across rules |
| RAM | 64 GB | ChExMix requests a 48 GB Java heap |
| Disk | 120 GB | inputs, intermediates, and outputs |
| Operating system | Linux or macOS | tested on RHEL 8 |
| Internet access | required for two rules only | see step 5 |

A laptop with 64 GB of RAM will work. On an HPC cluster, request an interactive or batch allocation with the resources above.

### 2. Install conda

If you do not already have conda or mamba, install Miniforge:

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
conda init bash
```

Close and reopen your terminal after `conda init`.

### 3. Clone the repository and create the Snakemake environment

```bash
git clone https://github.com/d1stadeyemi/perr_chipexo.git
cd perr_chipexo
conda env create -f workflow/envs/snakemake.yaml
conda activate snakemake
snakemake --version        # expect 9.19.0
```

This installs Snakemake only. Every bioinformatics tool is installed automatically, one conda environment per rule, the first time you run the workflow. You do not need to install any of them yourself.

### 4. Place the input data

Two files are required. Neither is included in this repository.

```bash
mkdir -p data
# 1. the sequencing core delivery, unmodified
cp /path/to/844_LBurcham_UTenn.zip data/
# 2. the RNA-seq differential expression workbook
cp /path/to/DEG_supplementary.xlsx data/
```

If your filenames differ, update the two matching lines in `config/config.yaml`:

```yaml
input_zip: "data/844_LBurcham_UTenn.zip"
integration:
  deg_excel: "data/DEG_supplementary.xlsx"
```

Nothing else in the configuration needs to change to reproduce the published results.

### 5. Download references and motif databases

These two rules need internet access. On clusters where compute nodes are offline, run them on a login or transfer node first. They take a few minutes.

```bash
snakemake --sdm conda --conda-frontend conda -c1 get_reference_all get_motif_databases
```

This downloads the two genome assemblies from NCBI and the MEME Suite prokaryotic motif collections into `resources/`.

### 6. Run the workflow

```bash
snakemake --sdm conda --conda-frontend conda -c8 all 2>&1 | tee run.log
```

The first execution also builds the 14 conda environments, which adds roughly 30 minutes. Subsequent runs reuse them.

To confirm what will run before committing to it, add `-n` for a dry run:

```bash
snakemake -n --sdm conda --conda-frontend conda -c8 all
```

### 7. Running on a SLURM cluster

A profile is provided for the ISAAC-NG cluster at the University of Tennessee. Copy it and edit the account, partition and QOS for your site:

```bash
cp -r profiles/isaac profiles/mysite
$EDITOR profiles/mysite/config.yaml
snakemake --workflow-profile profiles/mysite all
```

The profile submits each rule as a separate SLURM job, using the `threads` and `resources` declared in that rule. To run the whole workflow inside a single allocation instead, request a node with the resources in step 1 and use the plain command from step 6.

### 8. Verify the run

```bash
cat results/02_report/qc_scorecard.md
cat results/11_sensitivity/SENSITIVITY_REPORT.md
column -t -s $'\t' results/07_figures/regulon_gene_counts.tsv
```

Every number reported in the publication is present in `results/`. The tables under `results/06_integration/` and `results/08_orthology/` are the primary outputs.

### 9. Running individual stages

Every stage has a named target and can be run alone. Snakemake resolves upstream dependencies automatically.

```bash
snakemake --sdm conda -c8 qc            # quality control only
snakemake --sdm conda -c8 peaks         # through peak calling
snakemake --sdm conda -c8 motif         # through motif discovery
snakemake --sdm conda -c8 integration   # through RNA-seq integration
```

### 10. Reproducing the exact software environment

The steps above pin the primary tools through `workflow/envs/*.yaml`, which is sufficient for most purposes. Conda still resolves transitive dependencies at install time, so two installations months apart can differ in libraries that are never named explicitly.

For byte-identical environments, `workflow/envs/locks/` contains fully resolved package lists captured from the run that produced the published results. Each file records every package, version and build, including transitive dependencies.

```bash
# rebuild any single environment exactly as it was
conda create --name perr_meme --file workflow/envs/locks/meme.lock.txt
conda activate perr_meme
meme -version                # 5.5.9
```

To make Snakemake use these instead of resolving fresh, replace the contents of the corresponding `workflow/envs/*.yaml` before the first run.

Two caveats. The lockfiles are platform specific and were generated on `linux-64`; they will not resolve on macOS or ARM. And they record package URLs, so they depend on those builds remaining available from the conda channels.

Regenerate them after any environment change:

```bash
for f in .snakemake/conda/*.yaml; do
  h=$(basename "$f" .yaml); d=".snakemake/conda/$h"
  [ -d "$d/conda-meta" ] || continue
  for src in workflow/envs/*.yaml; do
    diff -q <(grep -v '^\s*#' "$f" | grep -v '^\s*$') \
            <(grep -v '^\s*#' "$src" | grep -v '^\s*$') >/dev/null 2>&1 \
      && conda list --prefix "$d" --explicit \
           > "workflow/envs/locks/$(basename $src .yaml).lock.txt" && break
  done
done
```

---

## Pipeline

| Stage | Module | Target | Primary tools |
|---|---|---|---|
| Extraction and sample discovery | `setup.smk` | `stage` | Python |
| Read and library quality control | `qc.smk` | `qc` | samtools, FastQC, ChIPexoQual, deepTools |
| Aggregated QC report | `report.smk` | `report` | MultiQC |
| Reference acquisition | `reference.smk` | `get_reference_all` | NCBI datasets CLI |
| Peak calling | `peaks.smk` | `peaks` | ChExMix |
| Depth-matched control | `downsample.smk` | `downsample` | samtools, ChExMix |
| Motif discovery and scanning | `motif.smk` | `motif` | MEME-ChIP, STREME, TOMTOM, FIMO, bedtools |
| Motif database comparison | `motifdb.smk` | `motifdb` | TOMTOM |
| RNA-seq integration | `integration.smk` | `integration` | Python, pandas, openpyxl |
| Figures | `figures.smk` | `figures` | R, ggplot2, patchwork |
| Cross-serotype orthology | `orthology.smk` | `orthology` | BLAST+ |
| Orthologous-position signal test | `pileup.smk` | `pileup` | MUMmer4, deepTools, pyBigWig |
| Genome browser assets | `igv.smk` | `igv` | samtools, deepTools, UCSC tools |
| Parameter sensitivity | `sensitivity.smk` | `sensitivity_cheap` | Python |
| Cross-caller concordance | `concordance.smk` | `concordance` | Python |
| Motif-less peak analysis | `motifless.smk` | `motifless` | MEME-ChIP, SEA, TOMTOM |

---

## Configuration

All parameters live in `config/config.yaml`, which is commented throughout with the rationale for each value. The most consequential settings:

| Parameter | Value | Effect |
|---|---|---|
| `peaks.chexmix_mappability` | 0.98 | fraction of the genome uniquely alignable, used in the background model |
| `peaks.chexmix_qvalue` | 0.05 | event significance threshold |
| `peaks.chexmix_minfold` | 1.5 | minimum enrichment over the scaled control |
| `motif.meme_maxw` | 30 | maximum motif width in de novo discovery |
| `motif.fimo_thresh` | 1e-4 | p-value threshold for genome-wide motif scanning |
| `integration.max_upstream` | 300 | promoter window for peak to gene assignment |
| `integration.operon_max_gap` | 50 | same-strand gap treated as co-transcribed |
| `integration.min_reps` | 2 | minimum replicate support for a peak to be assigned to a gene |
| `integration.deg_padj` | 0.05 | differential expression significance |
| `integration.deg_log2fc` | 0.585 | differential expression effect size, equivalent to 1.5 fold |

The `sensitivity` module measures how much each of these matters. It runs the affected part of the pipeline across a grid of alternatives, writing to `results/11_sensitivity/` without ever modifying `config/config.yaml`, and reports which outputs are stable and which depend on the threshold.

---

## Implementation notes

Several stages required adaptation because the tools involved carry defaults calibrated for mammalian genomes. These are documented here because they are the parts of the workflow most likely to be relevant to anyone applying similar tools to a small, densely covered genome.

**Mappability must match the organism.** ChExMix defaults to a mappability of 0.8, which reflects repeat content in a 3 Gb genome. A 2.1 Mb bacterial chromosome is almost entirely uniquely alignable, and the correct value is close to 0.98. This changes the background model and materially changes the peak set.

**ChIPexoQual requires subsampling on small genomes.** At full depth on a 2 Mb genome, read coverage is effectively continuous and the entire chromosome collapses into a single read island, leaving no regression to fit. The workflow subsamples every sample to a fixed 250,000 reads before island formation. This also makes the resulting metrics comparable across samples of differing depth.

**One motif matrix is applied to both genomes.** The two serotypes yield discoveries that differ in orientation, motif identifier, and rank. Scanning each genome with its own matrix would compare counts obtained with different matrices, so the strongest single discovery is extracted and applied symmetrically.

**Peak to gene assignment handles divergent promoters explicitly.** A peak is assigned to every gene whose transcription start site lies within the promoter window on either strand, then extended through same-strand operon members. Nearest-gene assignment would discard one of the two genes at every divergent promoter, and bacterial genomes are dense with them.

**Conservation is computed on two independent axes.** One uses binding alone, with no reference to expression; the other combines binding with expression response. Collapsing these into a single label conflates cistrome conservation with regulatory conservation, and the two answer different questions.

**Positional motif enrichment is computed in the workflow rather than by CentriMo.** CentriMo, SpaMo, and the HTML writer used by SEA all fail in the MEME 5.5.9 conda build with a template error. The underlying analyses complete and their tabular outputs are written first, so results remain recoverable. The workflow computes the motif offset directly relative to the ChExMix crosslink position, which is also a more meaningful reference point than the sequence midpoint that CentriMo uses.

**Peaks with single-replicate support are used for discovery but not for assignment.** Motif discovery uses every peak, because more sequences means more statistical power. Gene assignment requires support from at least two replicates, so downstream claims about which genes are bound rest on reproducible events. The threshold is configurable and its effect is quantified by the sensitivity module.

**A genome build assertion halts the run on mismatch.** The workflow reads the `@SQ` header from every BAM and compares it against the accession declared in the configuration. A mismatch aborts execution rather than producing a warning.

---

## Tools and dependencies

Managed automatically by Snakemake through per-rule conda environments. No manual installation is required.

| Tool | Version | Purpose |
|---|---|---|
| Snakemake | 9 | workflow engine |
| samtools | 1.21 | BAM handling, subsampling, indexing |
| FastQC | 0.12 | raw read quality |
| MultiQC | 1.21 | aggregated QC reporting |
| ChIPexoQual | 1.28 | ChIP-exo specific library metrics |
| deepTools | 3.5 | coverage tracks, fingerprint, correlation |
| ChExMix | 0.52 | ChIP-exo peak calling |
| MEME Suite | 5.5.9 | MEME, STREME, TOMTOM, FIMO, SEA, iupac2meme |
| bedtools | 2.31 | interval operations, sequence extraction |
| BLAST+ | 2.15 | reciprocal best hit orthology |
| MUMmer4 | 4.0 | whole-genome alignment for coordinate lifting |
| NCBI datasets CLI | 16 | reference genome retrieval |
| pyBigWig | 0.3 | coverage extraction |
| GEM | 3.4 | independent peak caller, concordance check (installed manually, see below) |
| R | 4.3 | figures, via ggplot2, patchwork, ggrepel |
| Python | 3.12 | analysis scripts, via pandas and openpyxl |

Full version pinning is in `workflow/envs/*.yaml`; fully resolved lockfiles are in `workflow/envs/locks/`. Literature citations for each tool are in [`CITATIONS.md`](CITATIONS.md).

GEM is the one tool not available through conda. It is a Java archive distributed by its authors and is used only for the cross-caller concordance check, which is not part of the main pipeline. To install it:

```bash
mkdir -p resources/gem && cd resources/gem
wget https://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz
wget https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_ChIP-exo.txt
tar xzf gem.v3.4.tar.gz
```

It requires Java 8 or later.

---

## Limitations

Constraints of the data and the workflow that affect interpretation. These are stated in full in the publication.

- **Protocol asymmetry between strains.** The sequencing core prepared the CJB111 samples with half the standard bead and antibody slurry. This is documented in the core report, not inferred. It affects signal recovery and therefore every cross-strain comparison. Modules `04_downsample` and `09_pileup` exist to bound its effect.
- **Single growth condition.** Cells were harvested unstressed. The workflow cannot distinguish a site that is never bound from one bound only under peroxide stress.
- **Epitope validation is genetic, not biochemical.** The tagged construct was confirmed by sequencing. No western blot was performed.
- **Annotation drift.** Thirteen gene identifiers present in the RNA-seq workbook are absent from the current RefSeq annotation, indicating the expression data were quantified against a different release. This affects a small number of accessory-genome assignments.
- **Alignment is inherited.** The workflow starts from the aligned BAMs supplied by the core and does not re-align from FASTQ. The core's parameters are recorded from the BAM headers, but the alignment itself is not independently reproduced.
- **Motif databases are taxonomically skewed.** The curated collections used for comparison are dominated by *Escherichia coli* and *Bacillus subtilis*. A species-specific regulator absent from those collections cannot be detected by that comparison.
- **The second peak caller is not a pipeline stage.** GEM required parameters calibrated to this dataset's coverage depth, so encoding those values in a rule would misrepresent them as general. The run is preserved as a documented script; only the concordance calculation is automated.

---

## Generated documentation

These files are produced by the workflow and regenerate on every run, so they cannot drift from the results they describe.

| File | Contents |
|---|---|
| `results/02_report/qc_scorecard.md` | one row per sample with advisory quality flags |
| `results/11_sensitivity/SENSITIVITY_REPORT.md` | which reported numbers are threshold-dependent |
| `results/13_motifless/MOTIFLESS_REPORT.md` | analysis of peaks lacking the discovered motif |
| `results/10_igv/README.md` | how to open the genome browser sessions |

`docs/QC_METRICS.md` is maintained by hand and explains what each quality metric measures and why the thresholds were chosen.

---

## Publication

*Citation to be added on acceptance.*

The publication describes the biology, the results, and their interpretation. This repository contains the analysis that produced them.

Sequencing was performed by the Cornell BRC Epigenomics Facility (RRID:SCR_021287).

---

## License

See [`LICENSE`](LICENSE).