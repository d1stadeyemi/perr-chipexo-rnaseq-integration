# =============================================================================
# figures.smk — publication figures from the integration tables
# =============================================================================
# Tables are hard to read; these panels carry the same information visually.
# ggplot2 + patchwork; each figure is 1 row x 2 panels.
#
# IMPORTANT counting note: the regulon table has one row per (peak, gene) pair,
# so a gene bound by several peaks appears several times. Every count in these
# figures is de-duplicated to UNIQUE GENES using a class precedence
# (DIRECT_PERBOX > DIRECT > BOUND_NOT_DE > BOUND_NO_RNA > INDIRECT), so a gene
# that is DIRECT at one peak is never also counted as BOUND_NOT_DE at another.
# =============================================================================

FIGS = os.path.join(RESULTS, "07_figures")


rule regulon_figures:
    """Four-panel figure per strain + a cross-strain composition bar."""
    input:
        regulon=expand(
            os.path.join(INTEG, "{strain}", "{strain}.direct_regulon.tsv"),
            strain=STRAINS),
        rnaseq=expand(
            os.path.join(INTEG, "{strain}", "{strain}.rnaseq.tsv"),
            strain=STRAINS),
    output:
        fig1=expand(os.path.join(FIGS, "{strain}_fig1_composition_volcano.pdf"),
                    strain=STRAINS),
        fig2=expand(os.path.join(FIGS, "{strain}_fig2_binding_targets.pdf"),
                    strain=STRAINS),
        composition=os.path.join(FIGS, "regulon_composition.pdf"),
        counts=os.path.join(FIGS, "regulon_gene_counts.tsv"),
    params:
        strains=STRAINS,
        outdir=FIGS,
    log:
        os.path.join(RESULTS, "logs", "regulon_figures.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/plot_regulon.R"

rule motif_evidence_figure:
    """The sequence case for the GBS motif: logo, positional enrichment, the
    like-for-like comparison against the S. pyogenes box, and occupancy.

    None of the four panels depends on a per-site FIMO q-value. Those run
    0.6-1.0 for a genome-wide 11-mer scan and are not the evidence."""
    input:
        classification=expand(
            os.path.join(MOTIF, "{strain}", "{strain}.all.peak_classification.tsv"),
            strain=STRAINS),
        pwm=os.path.join(MOTIF, "gbs_denovo.pwm.meme"),
        fimo_gbs=expand(os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.tsv"),
                        strain=STRAINS),
        fimo_gas=expand(os.path.join(MOTIF, "{strain}", "fimo", "fimo.tsv"),
                        strain=STRAINS),
    output:
        pdf=os.path.join(FIGS, "motif_evidence.pdf"),
        png=os.path.join(FIGS, "motif_evidence.png"),
    params:
        strains=STRAINS,
        halfwindow=config["motif"]["motif_in_peak_halfwindow"],
        outdir=FIGS,
    log:
        os.path.join(LOG, "motif_evidence_figure.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/plot_motif_evidence.R"


rule peak_position_figure:
    """Peak position relative to the nearest TSS, as a continuous distribution
    rather than three threshold-defined categories, with the genome's coding
    fraction as the null expectation for intragenic peaks."""
    input:
        assign=expand(
            os.path.join(INTEG, "{strain}", "{strain}.peak_gene_assignments.tsv"),
            strain=STRAINS),
        gff=[os.path.join(REF, strain_assembly(s), f"{strain_assembly(s)}.gff")
             for s in STRAINS],
    output:
        pdf=os.path.join(FIGS, "peak_position.pdf"),
        png=os.path.join(FIGS, "peak_position.png"),
    params:
        strains=STRAINS,
        max_upstream=config["integration"]["max_upstream"],
        outdir=FIGS,
    log:
        os.path.join(LOG, "peak_position_figure.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/plot_peak_position.R"


rule conservation_euler_figure:
    """Area-proportional Euler for cross-serotype binding and regulon overlap.

    The figure carries its own caveat: a Venn is read at a glance, and the
    strain-specific region is confounded by the declared half-slurry prep."""
    input:
        cross=os.path.join(ORTH, "crossstrain_regulon.tsv"),
    output:
        pdf=os.path.join(FIGS, "conservation_euler.pdf"),
        png=os.path.join(FIGS, "conservation_euler.png"),
    params:
        strains=STRAINS,
        outdir=FIGS,
    log:
        os.path.join(LOG, "conservation_euler_figure.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/plot_conservation_euler.R"


rule metaprofile_figure:
    """Average strand-separated 5' end signal aligned on every crosslink.

    The canonical ChIP-exo figure: a locus panel can be cherry-picked, an
    average over the whole peak set cannot. Split by motif status so direct and
    motif-less peaks can be compared, and shown per-peak as a heatmap so the
    reader sees the distribution rather than only its mean."""
    input:
        fwd=expand(os.path.join(IGV, "bigwig", "full_depth", "{strain}",
                                "{strain}.ip.fwd.cpm.bw"), strain=STRAINS),
        rev=expand(os.path.join(IGV, "bigwig", "full_depth", "{strain}",
                                "{strain}.ip.rev.cpm.bw"), strain=STRAINS),
        ctrl=expand(os.path.join(IGV, "bigwig", "full_depth", "{strain}",
                                 "{strain}.ctrl.cpm.bw"), strain=STRAINS),
        events=expand(os.path.join(PEAKS_CHEXMIX, "{strain}", "{strain}.events"),
                      strain=STRAINS),
        classification=expand(
            os.path.join(MOTIF, "{strain}", "{strain}.all.peak_classification.tsv"),
            strain=STRAINS),
        repcodes=expand(
            os.path.join(PEAKS_CHEXMIX, "{strain}", "{strain}.replicationcodes.table"),
            strain=STRAINS),
    output:
        pdf=os.path.join(FIGS, "metaprofile.pdf"),
        png=os.path.join(FIGS, "metaprofile.png"),
    params:
        strains=STRAINS,
        flank=config["figures"]["metaprofile_flank"],
        outdir=FIGS,
        min_reps=config["figures"]["metaprofile_min_reps"],
    log:
        os.path.join(LOG, "metaprofile_figure.log"),
    conda:
        "../envs/metaprofile.yaml"
    script:
        "../scripts/plot_metaprofile.py"

rule figures:
    input:
        os.path.join(FIGS, "regulon_composition.pdf"),
        os.path.join(FIGS, "regulon_gene_counts.tsv"),
        os.path.join(FIGS, "orthology_fig0_binding_axis.pdf"),
        os.path.join(FIGS, "orthology_fig1_conservation.pdf"),
        os.path.join(FIGS, "orthology_fig2_core_concordance.pdf"),
        os.path.join(FIGS, "enrichment_distribution.pdf"),
        os.path.join(FIGS, "enrichment_summary.tsv"),
        os.path.join(FIGS, "motif_evidence.pdf"),
        os.path.join(FIGS, "peak_position.pdf"),
        os.path.join(FIGS, "conservation_euler.pdf"),
        os.path.join(FIGS, "metaprofile.pdf"),

rule orthology_figures:
    """Cross-serotype conservation figures (1 row x 2 panels each)."""
    input:
        crossstrain=os.path.join(ORTH, "crossstrain_regulon.tsv"),
        summary=os.path.join(ORTH, "crossstrain_summary.tsv"),
    output:
        fig0=os.path.join(FIGS, "orthology_fig0_binding_axis.pdf"),
        fig1=os.path.join(FIGS, "orthology_fig1_conservation.pdf"),
        fig2=os.path.join(FIGS, "orthology_fig2_core_concordance.pdf"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
        outdir=FIGS,
    log:
        os.path.join(RESULTS, "logs", "orthology_figures.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/plot_orthology.R"

rule enrichment_figures:
    """Enrichment distribution of called events, stratified by subtype and by
    replication code.

    Kept in the workflow rather than run by hand so it cannot go stale: it reads
    the peak calls directly, and every parameter change that moves the peak set
    regenerates it."""
    input:
        events=expand(
            os.path.join(PEAKS_CHEXMIX, "{strain}", "{strain}.events"),
            strain=STRAINS),
        repcodes=expand(
            os.path.join(PEAKS_CHEXMIX, "{strain}", "{strain}.replicationcodes.table"),
            strain=STRAINS),
    output:
        pdf=os.path.join(FIGS, "enrichment_distribution.pdf"),
        png=os.path.join(FIGS, "enrichment_distribution.png"),
        tsv=os.path.join(FIGS, "enrichment_summary.tsv"),
    params:
        strains=STRAINS,
        outdir=FIGS,
    log:
        os.path.join(LOG, "enrichment_figures.log"),
    conda:
        "../envs/figures.yaml"
    script:
        "../scripts/enrichment_distribution.R"