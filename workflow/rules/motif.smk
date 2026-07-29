# =============================================================================
# motif.smk — motif discovery, known-Per-box comparison, and peak classification
# =============================================================================
# Logic (important):
#   * MEME-ChIP performs UNBIASED de novo discovery of motifs over-represented
#     under the ChExMix PerR peaks. This is the actual result.
#   * The GAS Per box consensus TTANAATNATTNTAA (group A Streptococcus; NO
#     GBS-specific Per box is defined — that is part of this project's novelty)
#     is used ONLY as a REFERENCE for two comparisons against that discovery:
#       - TOMTOM: does the discovered motif resemble the GAS Per box?
#       - FIMO:   how often does the GAS Per box pattern occur, under peaks and
#                 genome-wide?
#   The GAS box never drives discovery; it is the yardstick.
#
# All tool invocations are grounded in the MEME Suite / bedtools documentation.
# An AT-rich GBS Markov background is built once per assembly and reused by
# iupac2meme, meme-chip, and fimo for internal consistency.
# =============================================================================

MOTIF = os.path.join(RESULTS, "05_motif")

# Peak sets analysed in parallel:
#   all    = every ChExMix binding event
#   hiconf = events reproducible across >= motif.hiconf_min_reps replicates
PEAKSETS = ["all", "hiconf"]

wildcard_constraints:
    peakset="all|hiconf",


# ---------------------------------------------------------------------------
# 1. AT-rich background model (per assembly), reused everywhere
#    fasta-get-markov: -m sets Markov order; -dna sets alphabet; both strands
#    by default. Low order captures base composition so the genome's baseline
#    AT-richness is not mistaken for a motif.
# ---------------------------------------------------------------------------
rule motif_background:
    input:
        fna=os.path.join(REF, "{asm}", "{asm}.fna"),
    output:
        bg=os.path.join(MOTIF, "background", "{asm}.m{order}.bg"),
    params:
        order=lambda wc: wc.order,
    log:
        os.path.join(RESULTS, "logs", "motif_bg_{asm}_m{order}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        "fasta-get-markov -m {params.order} -dna {input.fna} > {output.bg} 2> {log}"


# ---------------------------------------------------------------------------
# 2. GAS Per box as a MEME-format motif (built once, background-corrected)
#    iupac2meme: direct IUPAC input, -bg background, -dna alphabet.
# ---------------------------------------------------------------------------
rule perbox_motif:
    input:
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
    output:
        meme=os.path.join(MOTIF, "{strain}", "GAS_perbox.reference.meme"),
    params:
        consensus=config["motif"]["gas_perbox"],
    log:
        os.path.join(RESULTS, "logs", "perbox_motif_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        "iupac2meme -dna -bg {input.bg} {params.consensus} > {output.meme} 2> {log}"

# ---------------------------------------------------------------------------
# 2b. The de novo GBS PerR motif as a MEME-format motif.
#     This is the ACTUAL result of the discovery step, re-expressed as a
#     scannable IUPAC consensus so both genomes are scanned with an IDENTICAL
#     matrix. Scanning each strain with its own PWM would compare counts
#     obtained with two different matrices.
# ---------------------------------------------------------------------------
rule gbs_perbox_motif:
    input:
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
    output:
        meme=os.path.join(MOTIF, "{strain}", "GBS_perbox.denovo.meme"),
    params:
        consensus=config["motif"]["gbs_perbox"],
    log:
        os.path.join(RESULTS, "logs", "gbs_perbox_motif_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        "iupac2meme -dna -bg {input.bg} {params.consensus} > {output.meme} 2> {log}"


# ---------------------------------------------------------------------------
# 7b. Genome-wide scan for the DE NOVO GBS motif (per assembly).
#     Same threshold and same background as the GAS scan, so the two hit counts
#     are directly comparable.
# ---------------------------------------------------------------------------
rule fimo_scan_gbs:
    input:
        motif=(os.path.join(MOTIF, "gbs_denovo.pwm.meme")
               if config["motif"]["gbs_motif_source"] == "pwm"
               else os.path.join(MOTIF, "{strain}", "GBS_perbox.denovo.meme")),
        fna=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
    output:
        tsv=os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.tsv"),
        gff=os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.gff"),
    params:
        outdir=os.path.join(MOTIF, "{strain}", "fimo_gbs"),
        thresh=config["motif"]["fimo_thresh"],
    log:
        os.path.join(RESULTS, "logs", "fimo_gbs_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        fimo --oc {params.outdir} --bfile {input.bg} \
            --thresh {params.thresh} \
            {input.motif} {input.fna} > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# 3. Fixed-width windows around each ChExMix event Point
#    MEME-ChIP wants equal-length, summit-centred sequences. We take
#    ±(window/2) around each event Point, clipped to chromosome ends.
# ---------------------------------------------------------------------------
rule events_to_bed:
    input:
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        repcodes=os.path.join(PEAKS, "{strain}", "{strain}.replicationcodes.table"),
        fai=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna.fai"),
    output:
        bed=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.windows.bed"),
    params:
        half=lambda wc: config["motif"]["window"] // 2,
        peakset=lambda wc: wc.peakset,
        min_reps=config["motif"]["hiconf_min_reps"],
    log:
        os.path.join(RESULTS, "logs", "events_to_bed_{strain}_{peakset}.log"),
    conda:
        "../envs/meme.yaml"
    script:
        "../scripts/events_to_bed.py"


# ---------------------------------------------------------------------------
# 4. Extract sequences under the windows (strand-aware)
#    bedtools getfasta: -s respects strand.
# ---------------------------------------------------------------------------
rule extract_peak_fasta:
    input:
        bed=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.windows.bed"),
        fna=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna"),
    output:
        fa=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.peaks.fa"),
    log:
        os.path.join(RESULTS, "logs", "extract_fasta_{strain}_{peakset}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        "bedtools getfasta -fi {input.fna} -bed {input.bed} -s > {output.fa} 2> {log}"


# ---------------------------------------------------------------------------
# 5. De novo discovery — the actual result
#    meme-chip: -db runs its built-in TOMTOM vs the GAS Per box reference;
#    -bfile the AT-rich background; -meme-mod zoops (one site per window);
#    width bounds bracket the 15 bp Per box.
# ---------------------------------------------------------------------------
rule meme_chip:
    input:
        fa=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.peaks.fa"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
        perbox=os.path.join(MOTIF, "{strain}", "GAS_perbox.reference.meme"),
    output:
        html=os.path.join(MOTIF, "{strain}", "memechip_{peakset}", "meme-chip.html"),
        done=touch(os.path.join(MOTIF, "{strain}", ".memechip_{peakset}.done")),
    params:
        outdir=os.path.join(MOTIF, "{strain}", "memechip_{peakset}"),
        minw=config["motif"]["meme_minw"],
        maxw=config["motif"]["meme_maxw"],
        mod=config["motif"]["meme_mod"],
        extra=config["motif"].get("memechip_extra", ""),
    log:
        os.path.join(RESULTS, "logs", "memechip_{strain}_{peakset}.log"),
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
            -db {input.perbox} \
            -bfile {input.bg} \
            -meme-mod {params.mod} \
            -minw {params.minw} -maxw {params.maxw} \
            -dna -meme-p {threads} \
            {params.extra} \
            {input.fa} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# 6. Explicit discovered-vs-known comparison
#    tomtom <query motifs> <target motifs>. Query = MEME-ChIP's combined
#    discovered motifs; target = the GAS Per box reference.
# ---------------------------------------------------------------------------
rule tomtom_perbox:
    input:
        memechip_done=os.path.join(MOTIF, "{strain}", ".memechip_{peakset}.done"),
        perbox=os.path.join(MOTIF, "{strain}", "GAS_perbox.reference.meme"),
    output:
        tsv=os.path.join(MOTIF, "{strain}", "tomtom_{peakset}", "tomtom.tsv"),
    params:
        memechip_dir=os.path.join(MOTIF, "{strain}", "memechip_{peakset}"),
        outdir=os.path.join(MOTIF, "{strain}", "tomtom_{peakset}"),
        thresh=config["motif"]["tomtom_thresh"],
    log:
        os.path.join(RESULTS, "logs", "tomtom_{strain}_{peakset}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        # MEME-ChIP writes combined discovered motifs to combined.meme
        q="{params.memechip_dir}/combined.meme"
        if [ ! -f "$q" ]; then
            q=$(find {params.memechip_dir} -name '*.meme' -not -name 'GAS_*' | head -1)
        fi
        tomtom -oc {params.outdir} -thresh {params.thresh} \
            "$q" {input.perbox} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# 7. Genome-wide GAS Per-box scan (per assembly)
#    fimo --bfile background; default p-value threshold from config.
# ---------------------------------------------------------------------------
rule fimo_scan:
    input:
        perbox=os.path.join(MOTIF, "{strain}", "GAS_perbox.reference.meme"),
        fna=lambda wc: os.path.join(
            REF, strain_assembly(wc.strain),
            f"{strain_assembly(wc.strain)}.fna"),
        bg=lambda wc: os.path.join(
            MOTIF, "background",
            f"{strain_assembly(wc.strain)}.m{config['motif']['bg_order']}.bg"),
    output:
        tsv=os.path.join(MOTIF, "{strain}", "fimo", "fimo.tsv"),
        gff=os.path.join(MOTIF, "{strain}", "fimo", "fimo.gff"),
    params:
        outdir=os.path.join(MOTIF, "{strain}", "fimo"),
        thresh=config["motif"]["fimo_thresh"],
    log:
        os.path.join(RESULTS, "logs", "fimo_{strain}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        fimo --oc {params.outdir} --bfile {input.bg} \
            --thresh {params.thresh} \
            {input.perbox} {input.fna} > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# 8. Classify each peak: canonical (GAS Per box under the footprint) vs not
#    Intersect ChExMix event windows with FIMO hits; join subtype from events.
# ---------------------------------------------------------------------------
rule classify_peaks:
    input:
        windows=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.windows.bed"),
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        fimo=os.path.join(MOTIF, "{strain}", "fimo", "fimo.tsv"),
        fimo_gbs=os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.tsv"),
    output:
        tsv=os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.peak_classification.tsv"),
    log:
        os.path.join(RESULTS, "logs", "classify_{strain}_{peakset}.log"),
    conda:
        "../envs/meme.yaml"
    script:
        "../scripts/classify_peaks.py"

rule gbs_pwm_motif:
    """Extract the de novo GBS PerR PWM from MEME-ChIP output.

    ONE matrix is used for BOTH genomes. Scanning each strain with its own PWM
    would compare counts obtained with two different matrices -- and the two
    discoveries are not even in the same orientation (A909 reported the reverse
    complement, CJB111 the forward strand), so their motif IDs and ranks differ.
    Taking the strongest single discovery and applying it symmetrically is the
    only clean cross-strain comparison."""
    input:
        combined=os.path.join(
            MOTIF, config["motif"]["gbs_pwm_strain"], "memechip_all", "combined.meme"),
    output:
        meme=os.path.join(MOTIF, "gbs_denovo.pwm.meme"),
    params:
        mid=config["motif"]["gbs_pwm_id"],
    log:
        os.path.join(RESULTS, "logs", "gbs_pwm_motif.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        meme-get-motif -id {params.mid} {input.combined} > {output.meme} 2> {log}
        grep -c '^MOTIF' {output.meme} >> {log}
        """

rule motif:
    """Aggregate target: motif analysis for every strain x peak set."""
    input:
        expand(os.path.join(MOTIF, "{strain}", "memechip_{peakset}", "meme-chip.html"),
               strain=STRAINS, peakset=PEAKSETS),
        expand(os.path.join(MOTIF, "{strain}", "tomtom_{peakset}", "tomtom.tsv"),
               strain=STRAINS, peakset=PEAKSETS),
        expand(os.path.join(MOTIF, "{strain}", "{strain}.{peakset}.peak_classification.tsv"),
               strain=STRAINS, peakset=PEAKSETS),
        expand(os.path.join(MOTIF, "{strain}", "fimo_gbs", "fimo.tsv"),
               strain=STRAINS),
