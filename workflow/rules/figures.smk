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


rule figures:
    """Aggregate target."""
    input:
        os.path.join(FIGS, "regulon_composition.pdf"),
        os.path.join(FIGS, "regulon_gene_counts.tsv"),
