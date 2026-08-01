# =============================================================================
# igv.smk — generate every asset needed to browse this dataset in IGV
# =============================================================================
# WHY THIS MODULE EXISTS
#   Genome-browser panels are evidence, not decoration, and they were the last
#   part of the analysis still being produced by hand. Hand-built panels are not
#   reproducible: track choice, data range, colour and window function all
#   affect what a reader concludes, and none of it is recorded anywhere.
#
#   This module makes them a pipeline product. One `snakemake igv` run emits
#   coverage tracks, annotation tracks, a locus lookup table, ready-to-open IGV
#   SESSION FILES with every display setting already applied, batch scripts for
#   headless snapshots, and a README documenting all of it.
#
# TWO DEPTHS, AND WHICH TO USE
#   full_depth     DEFAULT. Every read. Use for browsing, and for any figure
#                  about ONE strain. Depth-matching would discard ~57% of A909's
#                  reads and 23 of its 87 events for no benefit here.
#   depth_matched  Built from results/04_downsample/. Use ONLY when a figure
#                  makes an explicit A909-vs-CJB111 claim, where unequal depth
#                  would otherwise confound the comparison.
#   Both are generated; the session files say which is which in the track names
#   so the two cannot be confused after the fact.
#
# A NORMALISATION TRAP THIS MODULE AVOIDS
#   Filtering to one strand and then asking deepTools for CPM computes the scale
#   factor from THAT SUBSET, so forward and reverse end up on different scales
#   and the footprint looks asymmetric when it is not. Instead we compute ONE
#   scale factor from the unfiltered merged BAM and pass it to both strands via
#   --scaleFactor. The tracks are then directly comparable by construction.
#
# WHY 5'-END-ONLY COVERAGE
#   In ChIP-exo the 5' end of R1 IS the exonuclease stop, i.e. the footprint
#   border. --Offset 1 uses only that base. Full read coverage would smear a
#   ~10 bp feature across the whole read length. Set igv.offset: 0 in config to
#   fall back to full read coverage.
# =============================================================================

IGV = os.path.join(RESULTS, "10_igv")
IGV_DEPTHS = config["igv"]["depths"]


wildcard_constraints:
    depth="full_depth|depth_matched",
    strand="fwd|rev",


def igv_ip_bams(strain, depth):
    ids = strain_ids(strain, "IP")
    if depth == "full_depth":
        return [raw_bam(s) for s in ids]
    return [os.path.join(DS, "bam", s + ".ds.bam") for s in ids]


def igv_ctrl_bam(strain, depth):
    c = strain_ids(strain, "control")[0]
    return raw_bam(c) if depth == "full_depth" else os.path.join(DS, "bam", c + ".ds.bam")


def igv_bam_idx(paths):
    return [p + ".bai" for p in paths]


# ---------------------------------------------------------------------------
# 1. Merge the IP replicates, per strain per depth
# ---------------------------------------------------------------------------
rule igv_merge_ip:
    """Pooled IP signal. Merging rather than showing three tracks: for browsing,
    the question is where signal is, and pooled reads give the best-powered
    view. Per-replicate tracks are still emitted if igv.per_replicate is set."""
    input:
        bams=lambda wc: igv_ip_bams(wc.strain, wc.depth),
        bais=lambda wc: igv_bam_idx(igv_ip_bams(wc.strain, wc.depth)),
    output:
        bam=temp(os.path.join(IGV, "bam", "{depth}", "{strain}.ip.merged.bam")),
        bai=temp(os.path.join(IGV, "bam", "{depth}", "{strain}.ip.merged.bam.bai")),
    log:
        os.path.join(LOG, "igv_merge_{strain}_{depth}.log"),
    threads: 4
    conda:
        "../envs/igv.yaml"
    shell:
        r"""
        set -euo pipefail
        samtools merge -f -@ {threads} {output.bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {output.bam}
        echo "primary mapped: $(samtools view -c -F 0x904 {output.bam})" >> {log}
        """


# ---------------------------------------------------------------------------
# 2. Strand-split coverage of the 5' ends
# ---------------------------------------------------------------------------
rule igv_bigwig_ip_stranded:
    """One track per strand. R1's 5' end carries the exonuclease stop, so only
    R1 is used and it is split by the strand R1 mapped to:
        forward  --samFlagInclude 64 --samFlagExclude 16   (R1, not reverse)
        reverse  --samFlagInclude 80                       (R1 AND reverse, 64+16)
    The scale factor comes from the UNFILTERED merged BAM so both strands share
    one scale (see the header note)."""
    input:
        bam=os.path.join(IGV, "bam", "{depth}", "{strain}.ip.merged.bam"),
        bai=os.path.join(IGV, "bam", "{depth}", "{strain}.ip.merged.bam.bai"),
    output:
        bw=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.{strand}.cpm.bw"),
    params:
        binsize=config["igv"]["binsize"],
        offset=config["igv"]["offset"],
        include=lambda wc: 64 if wc.strand == "fwd" else 80,
        exclude=lambda wc: "--samFlagExclude 16" if wc.strand == "fwd" else "",
        offset_flag=lambda wc: (f"--Offset {config['igv']['offset']}"
                                if int(config["igv"]["offset"]) > 0 else ""),
    log:
        os.path.join(LOG, "igv_bw_{strain}_{depth}_{strand}.log"),
    threads: 4
    conda:
        "../envs/igv.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bw})
        n=$(samtools view -c -F 0x904 {input.bam})
        sf=$(awk -v n="$n" 'BEGIN{{ printf "%.10f", 1000000.0/n }}')
        echo "total primary reads = $n   shared scaleFactor = $sf" > {log}
        bamCoverage -b {input.bam} -o {output.bw} \
            --binSize {params.binsize} {params.offset_flag} \
            --samFlagInclude {params.include} {params.exclude} \
            --scaleFactor "$sf" --normalizeUsing None \
            -p {threads} >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# 3. Control track, same scale convention, not strand-split
# ---------------------------------------------------------------------------
rule igv_bigwig_control:
    """The no-tag control. Load this in EVERY panel: a peak shown without a flat
    control beneath it is not evidence, and it is the first thing a sceptical
    reader looks for."""
    input:
        bam=lambda wc: igv_ctrl_bam(wc.strain, wc.depth),
        bai=lambda wc: igv_ctrl_bam(wc.strain, wc.depth) + ".bai",
    output:
        bw=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ctrl.cpm.bw"),
    params:
        binsize=config["igv"]["binsize"],
        offset_flag=lambda wc: (f"--Offset {config['igv']['offset']}"
                                if int(config["igv"]["offset"]) > 0 else ""),
    log:
        os.path.join(LOG, "igv_bw_ctrl_{strain}_{depth}.log"),
    threads: 4
    conda:
        "../envs/igv.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bw})
        n=$(samtools view -c -F 0x904 {input.bam})
        sf=$(awk -v n="$n" 'BEGIN{{ printf "%.10f", 1000000.0/n }}')
        echo "total primary reads = $n   scaleFactor = $sf" > {log}
        bamCoverage -b {input.bam} -o {output.bw} \
            --binSize {params.binsize} {params.offset_flag} \
            --scaleFactor "$sf" --normalizeUsing None \
            -p {threads} >> {log} 2>&1
        """


# ---------------------------------------------------------------------------
# 4. Negated reverse strand, for the mirrored display
# ---------------------------------------------------------------------------
rule igv_bigwig_reverse_negated:
    """IGV cannot invert a bigWig, so we write a negated copy. Displayed beneath
    the forward track with range [-max, 0] it produces the canonical mirrored
    ChIP-exo footprint, which reads far more clearly than two stacked tracks."""
    input:
        bw=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.rev.cpm.bw"),
        fai=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.fna.fai"),
    output:
        bw=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.rev.neg.cpm.bw"),
    log:
        os.path.join(LOG, "igv_bw_neg_{strain}_{depth}.log"),
    conda:
        "../envs/igv.yaml"
    shell:
        r"""
        set -euo pipefail
        tmp=$(mktemp -d)
        cut -f1,2 {input.fai} > "$tmp/chrom.sizes"
        bigWigToBedGraph {input.bw} "$tmp/rev.bg" > {log} 2>&1
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,-$4}}' "$tmp/rev.bg" > "$tmp/neg.bg"
        bedGraphToBigWig "$tmp/neg.bg" "$tmp/chrom.sizes" {output.bw} >> {log} 2>&1
        rm -rf "$tmp"
        """


# ---------------------------------------------------------------------------
# 5. Annotation tracks — prefix-stripped and enriched
# ---------------------------------------------------------------------------
rule igv_tracks:
    """ChExMix writes 'chrNC_007432.1'; the FASTA, GFF and bigWigs use
    'NC_007432.1'. A track with the wrong prefix loads SILENTLY and displays
    nothing — it does not error. This rule normalises every annotation track and
    additionally writes peak BEDs carrying fold-enrichment, subtype and motif
    status in the name field, so hovering a peak in IGV shows its statistics."""
    input:
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        events_ds=os.path.join(DS, "peaks", "{strain}_ds", "{strain}_ds.events"),
        classif=os.path.join(MOTIF, "{strain}", "{strain}.all.peak_classification.tsv"),
        fimo_gbs=os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.gff"),
        fimo_gas=os.path.join(MOTIF, "{strain}", "fimo", "fimo.gff"),
        regulon=os.path.join(INTEG, "{strain}", "{strain}.direct_regulon.tsv"),
        gff=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.gff"),
    output:
        peaks=os.path.join(IGV, "tracks", "{strain}.peaks.bed"),
        peaks_ds=os.path.join(IGV, "tracks", "{strain}.peaks_depthmatched.bed"),
        motif_gbs=os.path.join(IGV, "tracks", "{strain}.motif_gbs.gff"),
        motif_gas=os.path.join(IGV, "tracks", "{strain}.motif_gas.gff"),
        regulon=os.path.join(IGV, "tracks", "{strain}.direct_regulon.bed"),
    params:
        strain=lambda wc: wc.strain,
    log:
        os.path.join(LOG, "igv_tracks_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/make_igv_tracks.py"


# ---------------------------------------------------------------------------
# 6. Locus lookup table — "I found an interesting gene, where do I look?"
# ---------------------------------------------------------------------------
rule igv_loci:
    """One row per gene touched by the analysis, with a ready-to-paste IGV
    coordinate string, the peaks assigned to it, motif status and expression.
    grep a locus tag or gene name and paste the coordinate into IGV's box."""
    input:
        regulon=os.path.join(INTEG, "{strain}", "{strain}.direct_regulon.tsv"),
        assign=os.path.join(INTEG, "{strain}", "{strain}.peak_gene_assignments.tsv"),
        gff=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain), f"{strain_assembly(wc.strain)}.gff"),
    output:
        tsv=os.path.join(IGV, "loci", "{strain}.loci.tsv"),
    params:
        strain=lambda wc: wc.strain,
        flank_context=config["igv"]["flank_context"],
        flank_zoom=config["igv"]["flank_zoom"],
    log:
        os.path.join(LOG, "igv_loci_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/make_igv_loci.py"


rule igv_flagship_loci:
    """The curated panel list: the loci that carry the paper, plus the artefact
    candidate and the sensitivity-test pair. Derived from the analysis outputs
    rather than hard-coded, so it stays correct if the results change."""
    input:
        loci=expand(os.path.join(IGV, "loci", "{strain}.loci.tsv"), strain=STRAINS),
        crossstrain=os.path.join(ORTH, "crossstrain_regulon.tsv"),
        orthoanchor=os.path.join(RESULTS, "09_pileup", "orthoanchor_signal.tsv"),
        events=expand(os.path.join(PEAKS, "{strain}", "{strain}.events"), strain=STRAINS),
    output:
        tsv=os.path.join(IGV, "loci", "flagship_loci.tsv"),
    params:
        strains=STRAINS,
        igvdir=IGV,
    log:
        os.path.join(LOG, "igv_flagship_loci.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/make_igv_flagship.py"

# ---------------------------------------------------------------------------
# 6b. Stage the reference into the IGV tree
# ---------------------------------------------------------------------------
rule igv_stage_genome:
    """Copy the reference FASTA, its index and the GFF into results/10_igv/.

    WHY COPY RATHER THAN REFERENCE IN PLACE: the session files carry paths
    relative to their own directory, so every resource must live under
    results/10_igv/ for the tree to be portable. With this rule,

        tar czf igv_bundle.tar.gz results/10_igv

    is a complete bundle that opens on any machine with nothing outside it to
    chase. Cost is ~2 MB FASTA + ~5 MB GFF per strain, which is nothing against
    the bigWigs.

    The .fai is copied too: IGV needs it and will otherwise try to write one
    next to the FASTA, which fails if the bundle is extracted read-only."""
    input:
        fna=lambda wc: os.path.join(REF, wc.asm, f"{wc.asm}.fna"),
        fai=lambda wc: os.path.join(REF, wc.asm, f"{wc.asm}.fna.fai"),
        gff=lambda wc: os.path.join(REF, wc.asm, f"{wc.asm}.gff"),
    output:
        fna=os.path.join(IGV, "genome", "{asm}", "{asm}.fna"),
        fai=os.path.join(IGV, "genome", "{asm}", "{asm}.fna.fai"),
        gff=os.path.join(IGV, "genome", "{asm}", "{asm}.gff"),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.fna})
        cp -f {input.fna} {output.fna}
        cp -f {input.fai} {output.fai}
        # Drop the whole-chromosome 'region' feature: it carries no Name
        # attribute, so IGV labels it ANONYMOUS and, spanning 2.1 Mb, that
        # label sits on top of every gene name in every view.
        # NOTE curly braces are format placeholders here, even inside comments.
        awk -F'\t' '$3 != "region" && $3 != "sequence_feature"' \
            {input.gff} > {output.gff}
        """

# ---------------------------------------------------------------------------
# 7. Session files, batch scripts and documentation
# ---------------------------------------------------------------------------
rule igv_session:
    """An IGV session is XML, so every display setting can be written by script:
    track order, colour, height, explicit data range, autoscale OFF, and window
    function = Maximum (IGV's default of Mean averages away the single-base
    spikes that ARE the ChIP-exo signal).

    Opening one file therefore reproduces a fully configured browser, and the
    settings are version-controlled rather than living in someone's memory."""
    input:
        fwd=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.fwd.cpm.bw"),
        rev=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.rev.cpm.bw"),
        neg=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ip.rev.neg.cpm.bw"),
        ctrl=os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ctrl.cpm.bw"),
        peaks=os.path.join(IGV, "tracks", "{strain}.peaks.bed"),
        peaks_ds=os.path.join(IGV, "tracks", "{strain}.peaks_depthmatched.bed"),
        motif_gbs=os.path.join(IGV, "tracks", "{strain}.motif_gbs.gff"),
        motif_gas=os.path.join(IGV, "tracks", "{strain}.motif_gas.gff"),
        regulon=os.path.join(IGV, "tracks", "{strain}.direct_regulon.bed"),
        flagship=os.path.join(IGV, "loci", "flagship_loci.tsv"),
        # staged INSIDE results/10_igv/ so the session's relative paths stay
        # inside the bundle
        fna=lambda wc: os.path.join(
            IGV, "genome", strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna"),
        fai=lambda wc: os.path.join(
            IGV, "genome", strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna.fai"),
        gff=lambda wc: os.path.join(
            IGV, "genome", strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.gff"),
    output:
        session=os.path.join(IGV, "session", "{strain}_{depth}.xml"),
        session_mirror=os.path.join(IGV, "session", "{strain}_{depth}_mirrored.xml"),
        batch=os.path.join(IGV, "batch", "{strain}_{depth}_panels.bat"),
    params:
        strain=lambda wc: wc.strain,
        depth=lambda wc: wc.depth,
        cfg=config["igv"],
        igvdir=IGV,
    log:
        os.path.join(LOG, "igv_session_{strain}_{depth}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/make_igv_session.py"


rule igv_readme:
    """Auto-generated documentation: every file produced, what it is, the exact
    load order, the display settings baked into the sessions, and the
    interpretation caveats. Regenerated on every run so it cannot drift."""
    input:
        sessions=expand(os.path.join(IGV, "session", "{strain}_{depth}.xml"),
                        strain=STRAINS, depth=IGV_DEPTHS),
        loci=expand(os.path.join(IGV, "loci", "{strain}.loci.tsv"), strain=STRAINS),
        flagship=os.path.join(IGV, "loci", "flagship_loci.tsv"),
    output:
        md=os.path.join(IGV, "README.md"),
    params:
        strains=STRAINS,
        depths=IGV_DEPTHS,
        cfg=config["igv"],
        igvdir=IGV,
        refdir=REF,
        assemblies={s: strain_assembly(s) for s in STRAINS},
    log:
        os.path.join(LOG, "igv_readme.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/make_igv_readme.py"


rule igv:
    """Aggregate target: everything needed to open this dataset in IGV."""
    input:
        expand(os.path.join(IGV, "bigwig", "{depth}", "{strain}",
                            "{strain}.ip.{strand}.cpm.bw"),
               depth=IGV_DEPTHS, strain=STRAINS, strand=["fwd", "rev"]),
        expand(os.path.join(IGV, "bigwig", "{depth}", "{strain}",
                            "{strain}.ip.rev.neg.cpm.bw"),
               depth=IGV_DEPTHS, strain=STRAINS),
        expand(os.path.join(IGV, "bigwig", "{depth}", "{strain}", "{strain}.ctrl.cpm.bw"),
               depth=IGV_DEPTHS, strain=STRAINS),
        expand(os.path.join(IGV, "session", "{strain}_{depth}.xml"),
               strain=STRAINS, depth=IGV_DEPTHS),
        expand(os.path.join(IGV, "batch", "{strain}_{depth}_panels.bat"),
               strain=STRAINS, depth=IGV_DEPTHS),
        os.path.join(IGV, "loci", "flagship_loci.tsv"),
        os.path.join(IGV, "README.md"),


rule igv_full_depth:
    """Convenience: only the full-depth assets, which is what browsing needs."""
    input:
        expand(os.path.join(IGV, "session", "{strain}_full_depth.xml"), strain=STRAINS),
        os.path.join(IGV, "README.md"),
