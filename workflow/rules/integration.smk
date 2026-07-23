# =============================================================================
# integration.smk — ChIP-exo x RNA-seq: the DIRECT PerR regulon
# =============================================================================
# The scientific payoff. RNA-seq alone cannot separate genes PerR sits on from
# genes that shift as a downstream consequence. Combining the two gives:
#
#   DIRECT / DIRECT_PERBOX  bound at the promoter AND differentially expressed
#   BOUND_NOT_DE            bound, expression unchanged
#   INDIRECT                differentially expressed, no binding -> downstream
#
# Peak->gene assignment deliberately does NOT use "nearest gene": the strongest
# site here (A909 874,861) lies in a divergent intergenic region, where a bound
# repressor occludes BOTH promoters. We assign a peak to every gene whose TSS is
# within `max_upstream` bp downstream of it, on either strand, then extend
# through co-transcribed operon members.
# =============================================================================

INTEG = os.path.join(RESULTS, "06_integration")


rule parse_degs:
    """Excel workbook -> tidy per-strain RNA-seq table (uses the RAW sheet)."""
    input:
        xlsx=config["integration"]["deg_excel"],
    output:
        tsv=os.path.join(INTEG, "{strain}", "{strain}.rnaseq.tsv"),
    params:
        strain=lambda wc: wc.strain,
        raw_sheet=lambda wc: config["integration"]["sheets"][wc.strain]["raw"],
        deg_sheet=lambda wc: config["integration"]["sheets"][wc.strain]["degs"],
        padj=config["integration"]["deg_padj"],
        log2fc=config["integration"]["deg_log2fc"],
    log:
        os.path.join(RESULTS, "logs", "parse_degs_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/parse_deg_excel.py"


rule assign_peaks:
    """Assign each binding event to gene(s) / transcription unit(s)."""
    input:
        gff=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.gff"),
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        repcodes=os.path.join(PEAKS, "{strain}", "{strain}.replicationcodes.table"),
        classification=os.path.join(
            MOTIF, "{strain}", "{strain}.all.peak_classification.tsv"),
    output:
        tsv=os.path.join(INTEG, "{strain}", "{strain}.peak_gene_assignments.tsv"),
    params:
        strain=lambda wc: wc.strain,
        max_upstream=config["integration"]["max_upstream"],
        operon_max_gap=config["integration"]["operon_max_gap"],
    log:
        os.path.join(RESULTS, "logs", "assign_peaks_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/assign_peaks_to_genes.py"


rule regulon_table:
    """Join binding with expression -> the direct regulon table."""
    input:
        assign=os.path.join(INTEG, "{strain}", "{strain}.peak_gene_assignments.tsv"),
        deg=os.path.join(INTEG, "{strain}", "{strain}.rnaseq.tsv"),
    output:
        tsv=os.path.join(INTEG, "{strain}", "{strain}.direct_regulon.tsv"),
        summary=os.path.join(INTEG, "{strain}", "{strain}.regulon_summary.tsv"),
    params:
        strain=lambda wc: wc.strain,
    log:
        os.path.join(RESULTS, "logs", "regulon_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/build_regulon_table.py"


rule integration:
    """Aggregate target."""
    input:
        expand(os.path.join(INTEG, "{strain}", "{strain}.direct_regulon.tsv"),
               strain=STRAINS),
        expand(os.path.join(INTEG, "{strain}", "{strain}.regulon_summary.tsv"),
               strain=STRAINS),
