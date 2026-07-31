# =============================================================================
# sensitivity.smk — how much does each threshold actually matter?
# =============================================================================
# WHY THIS MODULE EXISTS
#   Every headline number in this study depends on parameters that were chosen,
#   not measured: a 300 bp promoter window, a 50 bp operon gap, padj < 0.05,
#   |log2FC| >= 1, FIMO p < 1e-4, 30% identity for orthology, ChExMix q < 0.01
#   and --minfold 1.5. Each is defensible. None is inevitable.
#
#   A reviewer will ask how stable the answer is. "We did not check" is the only
#   answer that cannot be defended, and the check is cheap for most axes.
#
#   The first run of this analysis (done as a shell loop) already showed why it
#   matters: the per-strain direct-regulon count ranges 11-24 across the
#   assignment grid, while SHARED_DIRECT is invariant at 5. The conserved core
#   is a property of the data; the absolute count is a property of a threshold.
#   Only the sweep distinguishes them.
#
# THE DESIGN RULE THAT MATTERS
#   Everything is written to a PARALLEL TREE at results/11_sensitivity/.
#   config/config.yaml is NEVER modified.
#
#   The earlier shell loop mutated config.yaml in place and restored it
#   afterwards. The restore worked, but the OUTPUT FILES were left at the last
#   iteration's settings (400/100) while the config said 300/50. Snakemake
#   caught it on the next run, but only because the params trigger exists. This
#   module cannot create that state: the main pipeline's files are read, never
#   written.
#
# TIERS, BY COST
#   cheap      seconds       assignment, DEG thresholds, orthology floors
#   moderate   minutes       FIMO threshold (rescan, reclassify, reassign)
#   expensive  hours         MEME width/model, ChExMix q/minfold
#
#   Run `sensitivity_cheap` routinely. The expensive tiers are opt-in because
#   ChExMix takes ~28 min per strain per combination.
#
# SCOPE: WHAT EACH AXIS PROPAGATES
#   assignment, deg, orthology, fimo  -> full chain to crossstrain_summary, so
#                                        the effect on every headline number is
#                                        visible
#   meme, chexmix                     -> summarised at their own level only
#                                        (discovered motifs; peak counts).
#                                        Propagating a 12-combination ChExMix
#                                        grid through the whole pipeline would
#                                        cost days for a question those two
#                                        axes answer directly.
# =============================================================================

import itertools

SENS = os.path.join(RESULTS, "11_sensitivity")
_SC = config.get("sensitivity", {})


def _grid(axis, keys, fmt):
    """Expand one axis's parameter lists into {combo_name: {param: value}}.

    Values are kept as STRINGS in config and used verbatim in directory names,
    so that 1.0 and 1 cannot silently produce two different paths for the same
    setting."""
    spec = None
    for tier in ("cheap", "moderate", "expensive"):
        if axis in _SC.get(tier, {}):
            spec = _SC[tier][axis]
            break
    if not spec:
        return {}
    out = {}
    for combo in itertools.product(*[spec[k] for k in keys]):
        out[fmt.format(*combo)] = dict(zip(keys, combo))
    return out


SENS_GRID = {
    "assignment": _grid("assignment", ["max_upstream", "operon_max_gap"], "up{0}_gap{1}"),
    "deg":        _grid("deg", ["padj", "log2fc"], "padj{0}_lfc{1}"),
    "orthology":  _grid("orthology", ["min_identity", "min_coverage"], "pid{0}_cov{1}"),
    "fimo":       _grid("fimo", ["thresh"], "t{0}"),
    "meme":       _grid("meme", ["maxw", "mod"], "w{0}_{1}"),
    "chexmix":    _grid("chexmix", ["qvalue", "minfold"], "q{0}_mf{1}"),
}

PROPAGATING = ["assignment", "deg", "orthology", "fimo"]


def sens_combos(axis):
    return sorted(SENS_GRID.get(axis, {}))


def _p(axis, combo, key, default=None):
    """Parameter value for a combination, falling back to the main config so a
    rule shared between axes uses production settings on the axes that do not
    vary that parameter."""
    return SENS_GRID.get(axis, {}).get(combo, {}).get(key, default)


wildcard_constraints:
    axis="assignment|deg|orthology|fimo|meme|chexmix",
    combo="[A-Za-z0-9._+-]+",


# ---------------------------------------------------------------------------
# Input resolvers: use this axis's modified file where it exists, otherwise the
# main pipeline's. This is what lets one downstream rule serve every axis.
# ---------------------------------------------------------------------------
def sens_in_events(wc):
    if wc.axis == "chexmix":
        return os.path.join(SENS, "chexmix", wc.combo, f"{wc.strain}.events")
    return os.path.join(PEAKS, wc.strain, f"{wc.strain}.events")


def sens_in_classification(wc):
    if wc.axis == "fimo":
        return os.path.join(SENS, "fimo", wc.combo,
                            f"{wc.strain}.all.peak_classification.tsv")
    return os.path.join(MOTIF, wc.strain, f"{wc.strain}.all.peak_classification.tsv")


def sens_in_assign(wc):
    """peak_gene_assignments: regenerated for axes that change peak calling,
    peak classification or the assignment rule itself."""
    if wc.axis in ("assignment", "fimo", "chexmix"):
        return os.path.join(SENS, wc.axis, wc.combo,
                            f"{wc.strain}.peak_gene_assignments.tsv")
    return os.path.join(INTEG, wc.strain, f"{wc.strain}.peak_gene_assignments.tsv")


def sens_in_rnaseq(wc):
    if wc.axis == "deg":
        return os.path.join(SENS, "deg", wc.combo, f"{wc.strain}.rnaseq.tsv")
    return os.path.join(INTEG, wc.strain, f"{wc.strain}.rnaseq.tsv")


def sens_in_orthologs(wc):
    if wc.axis == "orthology":
        return os.path.join(SENS, "orthology", wc.combo, "orthologs_rbh.tsv")
    return os.path.join(ORTH, "orthologs_rbh.tsv")


# ===========================================================================
# CHEAP TIER
# ===========================================================================

rule sens_parse_degs:
    """DEG significance thresholds. The RAW sheet is re-read each time; only the
    padj and |log2FC| cutoffs change, so this measures how much of the direct /
    indirect split is a property of where the line was drawn."""
    input:
        xlsx=config["integration"]["deg_excel"],
    output:
        tsv=os.path.join(SENS, "deg", "{combo}", "{strain}.rnaseq.tsv"),
    params:
        strain=lambda wc: wc.strain,
        raw_sheet=lambda wc: config["integration"]["sheets"][wc.strain]["raw"],
        deg_sheet=lambda wc: config["integration"]["sheets"][wc.strain]["degs"],
        padj=lambda wc: _p("deg", wc.combo, "padj"),
        log2fc=lambda wc: _p("deg", wc.combo, "log2fc"),
    log:
        os.path.join(LOG, "sens_degs_{combo}_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/parse_deg_excel.py"


rule sens_ortholog_table:
    """Reciprocal-best-hit identity and coverage floors. BLAST is NOT re-run --
    the same blastp output is re-filtered, which is why this axis costs
    seconds rather than minutes."""
    input:
        fwd=os.path.join(ORTH, "blast",
                         f"{strain_assembly(STRAIN_A)}__vs__{strain_assembly(STRAIN_B)}.tsv"),
        rev=os.path.join(ORTH, "blast",
                         f"{strain_assembly(STRAIN_B)}__vs__{strain_assembly(STRAIN_A)}.tsv"),
        gff_a=os.path.join(REF, strain_assembly(STRAIN_A),
                           f"{strain_assembly(STRAIN_A)}.gff"),
        gff_b=os.path.join(REF, strain_assembly(STRAIN_B),
                           f"{strain_assembly(STRAIN_B)}.gff"),
    output:
        tsv=os.path.join(SENS, "orthology", "{combo}", "orthologs_rbh.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
        min_identity=lambda wc: _p("orthology", wc.combo, "min_identity"),
        min_coverage=lambda wc: _p("orthology", wc.combo, "min_coverage"),
        flag_identity=config["orthology"]["flag_identity"],
    log:
        os.path.join(LOG, "sens_orth_{combo}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/build_ortholog_table.py"


# ===========================================================================
# MODERATE TIER — FIMO threshold
# ===========================================================================

rule sens_fimo_gbs:
    """Genome-wide GBS-motif scan at a different p-value threshold. This is the
    axis that directly moves the '37 of 87 peaks carry the motif' figure, so it
    is worth knowing how much of that number is the threshold."""
    input:
        motif=(os.path.join(MOTIF, "gbs_denovo.pwm.meme")
               if config["motif"]["gbs_motif_source"] == "pwm"
               else os.path.join(MOTIF, "{strain}", "GBS_perbox.denovo.meme")),
        fna=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.fna"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
    output:
        tsv=os.path.join(SENS, "fimo", "{combo}", "{strain}", "fimo.tsv"),
        gff=os.path.join(SENS, "fimo", "{combo}", "{strain}", "fimo.gff"),
    params:
        outdir=os.path.join(SENS, "fimo", "{combo}", "{strain}"),
        thresh=lambda wc: _p("fimo", wc.combo, "thresh"),
    log:
        os.path.join(LOG, "sens_fimo_{combo}_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        fimo --oc {params.outdir} --bfile {input.bg} \
            --thresh {params.thresh} \
            {input.motif} {input.fna} > {log} 2>&1
        """


rule sens_classify_peaks:
    """Re-classify peaks against the rescanned motif hits. The GAS scan is left
    at production settings: only the GBS threshold is under test."""
    input:
        windows=os.path.join(MOTIF, "{strain}", "{strain}.all.windows.bed"),
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        fimo=os.path.join(MOTIF, "{strain}", "fimo", "fimo.tsv"),
        fimo_gbs=os.path.join(SENS, "fimo", "{combo}", "{strain}", "fimo.tsv"),
    output:
        tsv=os.path.join(SENS, "fimo", "{combo}",
                         "{strain}.all.peak_classification.tsv"),
    log:
        os.path.join(LOG, "sens_classify_{combo}_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    script:
        "../scripts/classify_peaks.py"


# ===========================================================================
# SHARED DOWNSTREAM — one chain serving every propagating axis
# ===========================================================================

rule sens_assign_peaks:
    input:
        gff=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.gff"),
        events=sens_in_events,
        repcodes=os.path.join(PEAKS, "{strain}", "{strain}.replicationcodes.table"),
        classification=sens_in_classification,
    output:
        tsv=os.path.join(SENS, "{axis}", "{combo}",
                         "{strain}.peak_gene_assignments.tsv"),
    params:
        strain=lambda wc: wc.strain,
        max_upstream=lambda wc: _p("assignment", wc.combo, "max_upstream",
                                   config["integration"]["max_upstream"]),
        operon_max_gap=lambda wc: _p("assignment", wc.combo, "operon_max_gap",
                                     config["integration"]["operon_max_gap"]),
    log:
        os.path.join(LOG, "sens_assign_{axis}_{combo}_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/assign_peaks_to_genes.py"


rule sens_regulon_table:
    input:
        assign=sens_in_assign,
        deg=sens_in_rnaseq,
        gff=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.gff"),
    output:
        tsv=os.path.join(SENS, "{axis}", "{combo}", "{strain}.direct_regulon.tsv"),
        summary=os.path.join(SENS, "{axis}", "{combo}", "{strain}.regulon_summary.tsv"),
    params:
        strain=lambda wc: wc.strain,
        perbox_source=config["integration"]["perbox_source"],
    log:
        os.path.join(LOG, "sens_regulon_{axis}_{combo}_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/build_regulon_table.py"


rule sens_crossstrain:
    input:
        reg_a=os.path.join(SENS, "{axis}", "{combo}", f"{STRAIN_A}.direct_regulon.tsv"),
        reg_b=os.path.join(SENS, "{axis}", "{combo}", f"{STRAIN_B}.direct_regulon.tsv"),
        orthologs=sens_in_orthologs,
    output:
        tsv=os.path.join(SENS, "{axis}", "{combo}", "crossstrain_regulon.tsv"),
        summary=os.path.join(SENS, "{axis}", "{combo}", "crossstrain_summary.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
    log:
        os.path.join(LOG, "sens_crossstrain_{axis}_{combo}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/crossstrain_regulon.py"


# ===========================================================================
# EXPENSIVE TIER
# ===========================================================================

rule sens_meme_chip:
    """Motif width and occurrence model.

    THE QUESTION: is TTMTAAATAAG a complete site, or one half of a longer
    inverted repeat? Fur-family proteins bind as dimers to a 7-1-7 core, and the
    production run capped width at 20 bp, which would truncate anything longer.

    -meme-mod also matters and answers a DIFFERENT question:
        zoops  zero or one occurrence per sequence (production setting)
        anr    any number of repetitions -- the only mode that can report TWO
               copies in one window, which is what a dimer on an inverted repeat
               would produce

    Caution: anr on AT-rich sequence over-reports, so read the E-values rather
    than the motif count. Summarised at motif level only -- propagating six
    combinations through the whole regulon would answer a question these two
    parameters do not ask."""
    input:
        fa=os.path.join(MOTIF, "{strain}", "{strain}.all.peaks.fa"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
        perbox=os.path.join(MOTIF, "{strain}", "GAS_perbox.reference.meme"),
    output:
        html=os.path.join(SENS, "meme", "{combo}", "{strain}", "meme-chip.html"),
        combined=os.path.join(SENS, "meme", "{combo}", "{strain}", "combined.meme"),
        summary=os.path.join(SENS, "meme", "{combo}", "{strain}", "summary.tsv"),
    params:
        outdir=os.path.join(SENS, "meme", "{combo}", "{strain}"),
        minw=config["motif"]["meme_minw"],
        maxw=lambda wc: _p("meme", wc.combo, "maxw"),
        mod=lambda wc: _p("meme", wc.combo, "mod"),
    log:
        os.path.join(LOG, "sens_memechip_{combo}_{strain}.log"),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240,
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        meme-chip -oc {params.outdir} \
            -db {input.perbox} -bfile {input.bg} \
            -meme-mod {params.mod} \
            -minw {params.minw} -maxw {params.maxw} \
            -dna -meme-p {threads} \
            {input.fa} > {log} 2>&1
        """


rule sens_chexmix:
    """Peak-calling significance and enrichment floors.

    --minfold is LINEAR fold, not log2: 1.0 disables the enrichment filter
    entirely, 1.5 is production, 2.0 is stricter. The weakest production A909
    event is log2FC 0.694 = 1.62-fold, just clearing the current floor, so the
    1.0 rung shows how many marginal events the filter currently removes.

    ~28 min per strain per combination. Opt-in for that reason."""
    input:
        bams=lambda wc: [raw_bam(s) for s in strain_ids(wc.strain)]
             + [raw_bai(s) for s in strain_ids(wc.strain)],
        geninfo=lambda wc: ref_geninfo(strain_assembly(wc.strain)),
    output:
        events=os.path.join(SENS, "chexmix", "{combo}", "{strain}.events"),
    params:
        outdir=os.path.join(SENS, "chexmix", "{combo}", "{strain}"),
        outname=lambda wc: wc.strain,
        signal=lambda wc: " ".join(
            f"--expt{wc.strain}-rep{i} {raw_bam(s)}"
            for i, s in enumerate(strain_ids(wc.strain, "IP"), start=1)),
        control=lambda wc: " ".join(
            f"--ctrl{wc.strain}-rep{i} {raw_bam(strain_ids(wc.strain,'control')[0])}"
            for i in range(1, len(strain_ids(wc.strain, "IP")) + 1)),
        fmt=config["peaks"]["format"],
        mappability=config["peaks"]["chexmix_mappability"],
        mem=32,
        qval=lambda wc: _p("chexmix", wc.combo, "qvalue"),
        minfold=lambda wc: _p("chexmix", wc.combo, "minfold"),
    log:
        os.path.join(LOG, "sens_chexmix_{combo}_{strain}.log"),
    threads: config["peaks"]["chexmix_threads"]
    resources:
        mem_mb=36000,
        runtime=600,
    conda:
        "../envs/chexmix.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        export _JAVA_OPTIONS="-Xmx{params.mem}G"
        chexmix --threads {threads} --geninfo {input.geninfo} \
            --format {params.fmt} {params.signal} {params.control} \
            --mappability {params.mappability} \
            --q {params.qval} --minfold {params.minfold} \
            --out {params.outdir}/{params.outname} > {log} 2>&1
        base="{params.outdir}/{params.outname}/{params.outname}_{params.outname}"
        if [ -f "${{base}}.events" ]; then cp "${{base}}.events" {output.events}
        else cp "$(find {params.outdir} -name '*.events' -not -path '*intermediate*' | head -1)" {output.events}
        fi
        """


# ===========================================================================
# SUMMARIES
# ===========================================================================

rule sens_summary_axis:
    """One row per combination, with every headline number recomputed."""
    wildcard_constraints:
        axis="assignment|deg|orthology|fimo",   # meme and chexmix have their own rules
    input:
        lambda wc: expand(
            os.path.join(SENS, wc.axis, "{combo}", "crossstrain_summary.tsv"),
            combo=sens_combos(wc.axis)),
        lambda wc: expand(
            os.path.join(SENS, wc.axis, "{combo}", "{strain}.direct_regulon.tsv"),
            combo=sens_combos(wc.axis), strain=STRAINS),
    output:
        tsv=os.path.join(SENS, "{axis}", "summary.tsv"),
    params:
        axis=lambda wc: wc.axis,
        combos=lambda wc: sens_combos(wc.axis),
        grid=lambda wc: SENS_GRID.get(wc.axis, {}),
        sensdir=SENS,
        strains=STRAINS,
    log:
        os.path.join(LOG, "sens_summary_{axis}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/sens_summary.py"


rule sens_summary_meme:
    input:
        expand(os.path.join(SENS, "meme", "{combo}", "{strain}", "summary.tsv"),
               combo=sens_combos("meme"), strain=STRAINS),
    output:
        tsv=os.path.join(SENS, "meme", "summary.tsv"),
    params:
        axis="meme",
        combos=sens_combos("meme"),
        grid=SENS_GRID.get("meme", {}),
        sensdir=SENS,
        strains=STRAINS,
    log:
        os.path.join(LOG, "sens_summary_meme.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/sens_summary.py"


rule sens_summary_chexmix:
    input:
        expand(os.path.join(SENS, "chexmix", "{combo}", "{strain}.events"),
               combo=sens_combos("chexmix"), strain=STRAINS),
    output:
        tsv=os.path.join(SENS, "chexmix", "summary.tsv"),
    params:
        axis="chexmix",
        combos=sens_combos("chexmix"),
        grid=SENS_GRID.get("chexmix", {}),
        sensdir=SENS,
        strains=STRAINS,
    log:
        os.path.join(LOG, "sens_summary_chexmix.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/sens_summary.py"


rule sens_report:
    """Markdown report: which numbers are robust, which are not, and the exact
    sentence to put in the manuscript for each."""
    input:
        lambda wc: [os.path.join(SENS, a, "summary.tsv")
                    for a in config["sensitivity"].get("report_axes",
                                                       ["assignment", "deg", "orthology"])],
    output:
        md=os.path.join(SENS, "SENSITIVITY_REPORT.md"),
    params:
        sensdir=SENS,
        axes=lambda wc: config["sensitivity"].get(
            "report_axes", ["assignment", "deg", "orthology"]),
        grids=SENS_GRID,
    log:
        os.path.join(LOG, "sens_report.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/sens_report.py"


# ===========================================================================
# TARGETS
# ===========================================================================

rule sensitivity_cheap:
    """Seconds. Run this routinely."""
    input:
        expand(os.path.join(SENS, "{axis}", "summary.tsv"),
               axis=[a for a in ("assignment", "deg", "orthology")
                     if sens_combos(a)]),


rule sensitivity_moderate:
    """Minutes. FIMO threshold."""
    input:
        os.path.join(SENS, "fimo", "summary.tsv"),


rule sensitivity_motif:
    """Hours. MEME width and occurrence model -- the half-site question."""
    input:
        os.path.join(SENS, "meme", "summary.tsv"),


rule sensitivity_peaks:
    """Many hours. ChExMix q-value and minfold. Opt-in."""
    input:
        os.path.join(SENS, "chexmix", "summary.tsv"),


rule sensitivity:
    """Everything except the ChExMix grid, plus the report."""
    input:
        rules.sensitivity_cheap.input,
        os.path.join(SENS, "fimo", "summary.tsv"),
        os.path.join(SENS, "SENSITIVITY_REPORT.md"),
