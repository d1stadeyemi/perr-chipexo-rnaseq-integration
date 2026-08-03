# =============================================================================
# pileup.smk — orthologous-position signal test
# =============================================================================
# THE QUESTION THIS MODULE EXISTS TO ANSWER
#   A909 has 85 bound genes, CJB111 has 28, and 25 are shared. CJB111's cistrome
#   is 89% NESTED inside A909's -- only 3 of its 28 bound genes are its own.
#   Depth-matching (results/04_downsample) removed sequencing depth as an
#   explanation: the ratio moved only from 3.8:1 to 3.2:1. Two candidates remain
#   and they are NOT distinguishable by any peak-calling statistic:
#
#     (i)  SENSITIVITY. CJB111 was prepared with half the standard bead/antibody
#          slurry (SignalFraction 0.102 vs 0.185). Lower enrichment means fewer
#          sites clear the significance threshold at ANY depth -- downsampling
#          scales signal and background together and cannot correct it.
#     (ii) BIOLOGY. The two serotypes genuinely have different PerR cistromes.
#
#   THE DISCRIMINATING TEST: go to the CJB111 position orthologous to each
#   A909-only site and ask what the read pileup looks like there, IGNORING
#   whether a peak was called.
#
#     elevated but sub-threshold  -> (i) sensitivity. The site is occupied;
#                                    the peak caller could not prove it.
#     flat background             -> (ii) biology. The site is genuinely unbound.
#
#   Both outcomes are publishable. The test is only meaningful against a NULL,
#   so signal is also extracted at orthologous positions of UNBOUND A909 genes;
#   that distribution is what "flat background" means quantitatively.
#
# TWO INDEPENDENT ANCHORING METHODS, KEPT IN SEPARATE COLUMNS AND NEVER MERGED
#   (a) ORTHOLOG-ANCHORED  Take the gene the A909 peak was assigned to, find its
#       CJB111 reciprocal best hit, and query a window at that gene's TSS.
#       Robust to indels and rearrangement; uses tables already built; fails for
#       events with no gene assignment or no ortholog.
#   (b) ALIGNMENT-LIFTED   nucmer-align the two chromosomes and lift the exact
#       A909 event coordinate into CJB111 space. Base-precise; handles events
#       with no gene assignment, including the divergent intergenic anchor at
#       A909 874,861; can fail across rearrangement breakpoints.
#   Where the two AGREE the result is strong. Where they DISAGREE that is
#   informative in itself -- it flags rearrangement or an assignment artifact --
#   which is exactly why they are not averaged into one number.
#
# COVERAGE IS BUILT FROM THE DEPTH-MATCHED BAMs (results/04_downsample/bam),
# NOT from the core's delivered bigwigs. The delivered bigwigs are full-depth
# with unknown normalisation; using them would silently reintroduce the depth
# confound that module 04 exists to remove.
# =============================================================================

PILE = os.path.join(RESULTS, "09_pileup")


def _ds_bam(sample):
    return os.path.join(DS, "bam", sample + ".ds.bam")


rule pileup_merge_ip:
    """One merged IP BAM per strain from the DEPTH-MATCHED replicates.

    Merging rather than averaging three separate bigwigs: at a candidate site
    the question is whether ANY consistent signal exists above background, and
    pooled reads give the best-powered single estimate of that. Per-replicate
    tracks remain available in results/04_downsample/bam if a reviewer asks for
    reproducibility across replicates."""
    input:
        bams=lambda wc: [_ds_bam(s) for s in strain_ids(wc.strain, "IP")],
        bais=lambda wc: [_ds_bam(s) + ".bai" for s in strain_ids(wc.strain, "IP")],
    output:
        bam=os.path.join(PILE, "bam", "{strain}.ip.merged.bam"),
        bai=os.path.join(PILE, "bam", "{strain}.ip.merged.bam.bai"),
    log:
        os.path.join(LOG, "pileup_merge_{strain}.log"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        samtools view -c -F 0x904 {output.bam} >> {log}
        """


rule pileup_bigwig_ip:
    """CPM-normalised, 1 bp resolution coverage of the merged IP.

    CPM (not RPKM/BPM): the two strains are already at matched read depth, so
    CPM keeps the tracks on a common scale without introducing a length term
    that is meaningless for point-source binding.

    binSize 1 because a ChIP-exo footprint is tens of bp wide; binning would
    smear the very feature under test."""
    input:
        bam=os.path.join(PILE, "bam", "{strain}.ip.merged.bam"),
        bai=os.path.join(PILE, "bam", "{strain}.ip.merged.bam.bai"),
    output:
        bw=os.path.join(PILE, "bigwig", "{strain}.ip.cpm.bw"),
    params:
        extra=config["pileup"]["bamcoverage_extra"],
    log:
        os.path.join(LOG, "pileup_bw_ip_{strain}.log"),
    threads: 4
    conda:
        "../envs/pileup.yaml"
    shell:
        r"""
        set -euo pipefail
        bamCoverage -b {input.bam} -o {output.bw} \
            --normalizeUsing CPM --binSize 1 \
            --numberOfProcessors {threads} {params.extra} > {log} 2>&1
        """


rule pileup_bigwig_ctrl:
    """Same treatment for the no-tag control, from its depth-matched BAM."""
    input:
        bam=lambda wc: _ds_bam(strain_ids(wc.strain, "control")[0]),
        bai=lambda wc: _ds_bam(strain_ids(wc.strain, "control")[0]) + ".bai",
    output:
        bw=os.path.join(PILE, "bigwig", "{strain}.ctrl.cpm.bw"),
    params:
        extra=config["pileup"]["bamcoverage_extra"],
    log:
        os.path.join(LOG, "pileup_bw_ctrl_{strain}.log"),
    threads: 4
    conda:
        "../envs/pileup.yaml"
    shell:
        r"""
        set -euo pipefail
        bamCoverage -b {input.bam} -o {output.bw} \
            --normalizeUsing CPM --binSize 1 \
            --numberOfProcessors {threads} {params.extra} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# METHOD (a): ortholog-anchored
# ---------------------------------------------------------------------------
rule pileup_orthoanchor:
    input:
        orthologs=os.path.join(ORTH, "orthologs_rbh.tsv"),
        assign_a=os.path.join(INTEG, STRAIN_A, f"{STRAIN_A}.peak_gene_assignments.tsv"),
        assign_b=os.path.join(INTEG, STRAIN_B, f"{STRAIN_B}.peak_gene_assignments.tsv"),
        gff_a=os.path.join(REF, strain_assembly(STRAIN_A),
                           f"{strain_assembly(STRAIN_A)}.gff"),
        gff_b=os.path.join(REF, strain_assembly(STRAIN_B),
                           f"{strain_assembly(STRAIN_B)}.gff"),
        bw_a_ip=os.path.join(PILE, "bigwig", f"{STRAIN_A}.ip.cpm.bw"),
        bw_a_ct=os.path.join(PILE, "bigwig", f"{STRAIN_A}.ctrl.cpm.bw"),
        bw_b_ip=os.path.join(PILE, "bigwig", f"{STRAIN_B}.ip.cpm.bw"),
        bw_b_ct=os.path.join(PILE, "bigwig", f"{STRAIN_B}.ctrl.cpm.bw"),
    output:
        tsv=os.path.join(PILE, "orthoanchor_signal.tsv"),
        summary=os.path.join(PILE, "orthoanchor_summary.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
        upstream=config["pileup"]["tss_upstream"],
        downstream=config["pileup"]["tss_downstream"],
    log:
        os.path.join(LOG, "pileup_orthoanchor.log"),
    conda:
        "../envs/pileup.yaml"
    script:
        "../scripts/pileup_orthoanchor.py"


# ---------------------------------------------------------------------------
# METHOD (b): whole-genome alignment lift
# ---------------------------------------------------------------------------
rule pileup_nucmer:
    """Whole-chromosome alignment, A909 (reference) vs CJB111 (query).

    --maxmatch keeps all anchor matches rather than only unique ones, which is
    the right choice for two strains of one species carrying repeated elements;
    delta-filter -1 then enforces a 1-to-1 mapping so each reference position
    lifts to at most one query position."""
    input:
        ref=os.path.join(REF, strain_assembly(STRAIN_A),
                         f"{strain_assembly(STRAIN_A)}.fna"),
        qry=os.path.join(REF, strain_assembly(STRAIN_B),
                         f"{strain_assembly(STRAIN_B)}.fna"),
    output:
        coords=os.path.join(PILE, "nucmer", "A_vs_B.1coords"),
    params:
        prefix=os.path.join(PILE, "nucmer", "A_vs_B"),
        minmatch=config["pileup"]["nucmer_minmatch"],
        mincluster=config["pileup"]["nucmer_mincluster"],
    log:
        os.path.join(LOG, "pileup_nucmer.log"),
    threads: 4
    conda:
        "../envs/pileup.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        nucmer --maxmatch -l {params.minmatch} -c {params.mincluster} \
            -t {threads} -p {params.prefix} {input.ref} {input.qry} > {log} 2>&1
        delta-filter -1 {params.prefix}.delta > {params.prefix}.1delta 2>> {log}
        show-coords -T -H -r {params.prefix}.1delta > {output.coords} 2>> {log}
        wc -l {output.coords} >> {log}
        """


rule pileup_lift:
    input:
        coords=os.path.join(PILE, "nucmer", "A_vs_B.1coords"),
        events_a=os.path.join(PEAKS_CHEXMIX, STRAIN_A, f"{STRAIN_A}.events"),
        events_b=os.path.join(PEAKS_CHEXMIX, STRAIN_B, f"{STRAIN_B}.events"),
        bw_a_ip=os.path.join(PILE, "bigwig", f"{STRAIN_A}.ip.cpm.bw"),
        bw_b_ip=os.path.join(PILE, "bigwig", f"{STRAIN_B}.ip.cpm.bw"),
        bw_b_ct=os.path.join(PILE, "bigwig", f"{STRAIN_B}.ctrl.cpm.bw"),
    output:
        tsv=os.path.join(PILE, "lifted_signal.tsv"),
        summary=os.path.join(PILE, "lifted_summary.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
        halfwin=config["pileup"]["lift_halfwindow"],
        n_null=config["pileup"]["n_null_positions"],
        seed=config["pileup"]["null_seed"],
    log:
        os.path.join(LOG, "pileup_lift.log"),
    conda:
        "../envs/pileup.yaml"
    script:
        "../scripts/pileup_lift.py"


rule pileup:
    """Aggregate target."""
    input:
        os.path.join(PILE, "orthoanchor_signal.tsv"),
        os.path.join(PILE, "orthoanchor_summary.tsv"),
        os.path.join(PILE, "lifted_signal.tsv"),
        os.path.join(PILE, "lifted_summary.tsv"),
