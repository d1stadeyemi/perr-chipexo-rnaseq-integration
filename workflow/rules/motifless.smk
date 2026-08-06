# =============================================================================
# motifless.smk — how does PerR recognise peaks that carry no GBS motif?
# =============================================================================
# THE QUESTION
#   Roughly 57% of A909 peaks (50 of 87) carry no GBS PerR motif at FIMO
#   p<1e-4. That is not anomalous -- almost every ChIP dataset contains peaks
#   lacking the canonical motif -- but the reason matters, and there are three
#   candidate explanations that make different predictions:
#
#     (a) DEGENERATE SITES. The motif is there but below threshold. Prediction:
#         de novo discovery on the motif-less subset recovers a variant of
#         TTMTAAATAAG.
#     (b) A SECOND BINDING MODE. PerR recognises a different sequence at these
#         sites. Prediction: de novo recovers a DIFFERENT motif.
#     (c) INDIRECT / TETHERED BINDING. PerR is brought to these sites by a
#         partner protein. Prediction: a partner TF's motif is over-represented,
#         and no PerR-like motif is found.
#
#   Only (a) and (b) are testable from sequence alone. (c) is supported, not
#   proven, by finding a partner motif -- Worsley Hunt & Wasserman (Genome Biol
#   2014) showed a small set of TF motifs recur as systematic artifacts across
#   hundreds of ChIP datasets, so any hit needs checking against that before it
#   is called a partner.
#
# WHY THE MOTIF-POSITIVE SET IS ALSO ANALYSED
#   A motif enriched in the motif-less peaks means nothing unless it is ALSO
#   checked in the motif-positive peaks. If it appears in both, it is a property
#   of PerR-bound regions generally (or of the peak-calling), not of the
#   motif-less subset. Both subsets therefore go through the identical pipeline
#   and the comparison is reported side by side.
#
# CAVEAT ON POWER
#   A909 has ~50 motif-less peaks and CJB111 ~8. MEME on 8 sequences is not
#   meaningful; the rule runs anyway but the summary flags any subset below
#   `min_seqs` as underpowered rather than reporting a motif from it.
# =============================================================================

ML = os.path.join(RESULTS, "13_motifless")


wildcard_constraints:
    subset="motifless|motifpos",


rule motifless_split:
    """Partition peak windows by GBS-motif status.

    The split uses the SAME classification table that drives DIRECT_PERBOX, so
    the subsets are defined identically to everything else in the analysis --
    no second definition of 'has a motif' to drift out of sync."""
    input:
        classif=os.path.join(MOTIF, "{strain}", "{strain}.all.peak_classification.tsv"),
        windows=os.path.join(MOTIF, "{strain}", "{strain}.all.windows.bed"),
    output:
        bed_neg=os.path.join(ML, "{strain}", "{strain}.motifless.bed"),
        bed_pos=os.path.join(ML, "{strain}", "{strain}.motifpos.bed"),
        tsv=os.path.join(ML, "{strain}", "{strain}.split_summary.tsv"),
    params:
        strain=lambda wc: wc.strain,
    log:
        os.path.join(LOG, "motifless_split_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/motifless_split.py"


rule motifless_fasta:
    input:
        bed=os.path.join(ML, "{strain}", "{strain}.{subset}.bed"),
        fna=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.fna"),
    output:
        fa=os.path.join(ML, "{strain}", "{strain}.{subset}.fa"),
    log:
        os.path.join(LOG, "motifless_fasta_{strain}_{subset}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        "bedtools getfasta -fi {input.fna} -bed {input.bed} -s > {output.fa} 2> {log}"


rule motifless_memechip:
    """4a: de novo discovery on each subset.

    -meme-nmotifs 5 rather than the default 1: the whole point is to see
    whether ANY sequence signal exists here, so several candidates are wanted.
    Background and width bounds are inherited from the main motif config so the
    comparison with the primary discovery is like for like."""
    input:
        fa=os.path.join(ML, "{strain}", "{strain}.{subset}.fa"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
        gbs=os.path.join(MOTIF, "gbs_denovo.pwm.meme"),
    output:
        html=os.path.join(ML, "{strain}", "memechip_{subset}", "meme-chip.html"),
        combined=os.path.join(ML, "{strain}", "memechip_{subset}", "combined.meme"),
        summary=os.path.join(ML, "{strain}", "memechip_{subset}", "summary.tsv"),
    params:
        outdir=os.path.join(ML, "{strain}", "memechip_{subset}"),
        minw=config["motif"]["meme_minw"],
        maxw=config["motif"]["meme_maxw"],
        mod=config["motif"]["meme_mod"],
        nmotifs=config["motifless"]["nmotifs"],
    log:
        os.path.join(LOG, "motifless_memechip_{strain}_{subset}.log"),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240,
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        # -db is the DISCOVERED GBS motif, so MEME-ChIP's own TOMTOM step
        # immediately answers "is what I found here just the GBS motif again?"
        meme-chip -oc {params.outdir} \
            -db {input.gbs} -bfile {input.bg} \
            -meme-mod {params.mod} -meme-nmotifs {params.nmotifs} \
            -minw {params.minw} -maxw {params.maxw} \
            -dna -meme-p {threads} \
            {input.fa} > {log} 2>&1
        """


rule motifless_tomtom_db:
    """4b: discovered motifs vs the 558 curated prokaryotic motifs.

    Same target database and threshold as motifdb.smk, so this result is
    directly comparable to the primary motif's negative result against the same
    reference set."""
    input:
        combined=os.path.join(ML, "{strain}", "memechip_{subset}", "combined.meme"),
        # config lists BARE FILENAMES; the directory comes from motifdb.smk
        dbs=[os.path.join(MDB, "PROKARYOTE", d)
             for d in config["motifdb"]["databases"]],
    output:
        tsv=os.path.join(ML, "{strain}", "tomtom_db_{subset}", "tomtom.tsv"),
    params:
        outdir=os.path.join(ML, "{strain}", "tomtom_db_{subset}"),
        thresh=config["motifdb"]["thresh"],
        dist=config["motifdb"]["dist"],
    log:
        os.path.join(LOG, "motifless_tomtomdb_{strain}_{subset}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        tomtom -oc {params.outdir} -thresh {params.thresh} \
            -dist {params.dist} {input.combined} {input.dbs} > {log} 2>&1
        """


rule motifless_report:
    """Side-by-side comparison and the interpretation each outcome supports."""
    input:
        splits=expand(os.path.join(ML, "{strain}", "{strain}.split_summary.tsv"),
                      strain=STRAINS),
        summaries=expand(
            os.path.join(ML, "{strain}", "memechip_{subset}", "summary.tsv"),
            strain=STRAINS, subset=["motifless", "motifpos"]),
        tomtom=expand(
            os.path.join(ML, "{strain}", "tomtom_db_{subset}", "tomtom.tsv"),
            strain=STRAINS, subset=["motifless", "motifpos"]),
    output:
        md=os.path.join(ML, "MOTIFLESS_REPORT.md"),
        tsv=os.path.join(ML, "motifless_summary.tsv"),
    params:
        strains=STRAINS,
        mldir=ML,
        min_seqs=config["motifless"]["min_seqs"],
    log:
        os.path.join(LOG, "motifless_report.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/motifless_report.py"


rule motifless:
    """Aggregate target."""
    input:
        os.path.join(ML, "MOTIFLESS_REPORT.md"),
