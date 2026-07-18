# =============================================================================
# peaks.smk — ChExMix peak (binding-event) calling
# =============================================================================
# ChExMix (Yamada et al. 2019; Pugh/Mahony) is the primary caller: a ChIP-exo
# mixture model, control-aware, reporting per-event signal/control tag counts
# and log2 fold. One run per strain: 3 PerR-His IP replicates vs that strain's
# no-tag control.
#
# PARAMETER CHOICES (all grounded in the ChExMix manual):
#   --mappability 0.98  Background/read-limit is Poisson-parameterised by
#                       reads / (mappable bases). Default 0.8 is a human value;
#                       for a 2.1 Mb genome with unique mapping ~0.98 is correct.
#                       This affects alpha and the per-base read limit, hence
#                       event counts — the single most important fix.
#   --format BAM        Documented token for BAM input (was SAM).
#   --q 0.01            Event significance (default; stated explicitly).
#   --minfold 1.5       Minimum fold over scaled control (default; explicit).
#   NCIS scaling        Default control normalisation (not overridden).
#   reporting: standard Default; conservative (significant in the condition as a
#                       whole). Kept deliberately; stated for methods.
#   --exclude           Optional artifact-region BED (rRNA/tRNA/prophage). Hook
#                       wired but empty until that BED rule is added; the manual
#                       recommends excluding blacklisted accumulations.
#
# Motif-defined subtypes (--seq/--back) are deferred to the motif stage, where
# an AT-rich GBS Markov background is built explicitly.
#
# Output layout (observed from the completed run): with --out {dir}/{name},
# ChExMix writes {dir}/{name}/{name}_{name}.events (+ .bed, replication table).
# =============================================================================


def chexmix_signal_flags(strain):
    """--exptSTRAIN-repN <bam> for each IP replicate."""
    return " ".join(
        f"--expt{strain}-rep{i} {raw_bam(s)}"
        for i, s in enumerate(strain_ids(strain, "IP"), start=1)
    )


def chexmix_control_flags(strain):
    """Map the single no-tag control to every signal replicate of the strain."""
    ctrl = strain_ids(strain, "control")
    if not ctrl:
        return ""
    ctrl_bam = raw_bam(ctrl[0])
    n_ip = len(strain_ids(strain, "IP"))
    return " ".join(f"--ctrl{strain}-rep{i} {ctrl_bam}" for i in range(1, n_ip + 1))


def chexmix_input_bams(strain):
    ids = strain_ids(strain)
    return [raw_bam(s) for s in ids] + [raw_bai(s) for s in ids]


rule chexmix:
    """Call PerR binding events for one strain (3 IP reps + no-tag control)."""
    input:
        bams=lambda wc: chexmix_input_bams(wc.strain),
        geninfo=lambda wc: ref_geninfo(strain_assembly(wc.strain)),
    output:
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        bed=os.path.join(PEAKS, "{strain}", "{strain}.bed"),
        repcodes=os.path.join(PEAKS, "{strain}", "{strain}.replicationcodes.table"),
        done=touch(os.path.join(PEAKS, "{strain}", ".chexmix.done")),
    params:
        outdir=os.path.join(PEAKS, "{strain}"),
        outname=lambda wc: wc.strain,
        signal=lambda wc: chexmix_signal_flags(wc.strain),
        control=lambda wc: chexmix_control_flags(wc.strain),
        fmt=config["peaks"]["format"],
        qval=config["peaks"]["chexmix_qvalue"],
        minfold=config["peaks"]["chexmix_minfold"],
        mappability=config["peaks"]["chexmix_mappability"],
        mem=config["peaks"]["chexmix_mem_gb"],
        exclude=lambda wc: (
            f"--exclude {config['peaks']['chexmix_exclude']}"
            if config["peaks"].get("chexmix_exclude") else ""
        ),
        extra=config["peaks"]["chexmix_extra"],
    log:
        os.path.join(RESULTS, "logs", "chexmix_{strain}.log"),
    threads: config["peaks"]["chexmix_threads"]
    resources:
        mem_mb=lambda wc: (config["peaks"]["chexmix_mem_gb"] + 4) * 1000,
        runtime=600,
    conda:
        "../envs/chexmix.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        export _JAVA_OPTIONS="-Xmx{params.mem}G"
        chexmix \
            --threads {threads} \
            --geninfo {input.geninfo} \
            --format {params.fmt} \
            {params.signal} \
            {params.control} \
            --mappability {params.mappability} \
            --q {params.qval} \
            --minfold {params.minfold} \
            {params.exclude} \
            --out {params.outdir}/{params.outname} \
            {params.extra} > {log} 2>&1

        # ChExMix writes into {outdir}/{name}/ with files prefixed {name}_{name}.
        base="{params.outdir}/{params.outname}/{params.outname}_{params.outname}"

        # events (required)
        if [ -f "${{base}}.events" ]; then
            cp "${{base}}.events" {output.events}
        else
            ev=$(find {params.outdir} -name '*.events' -not -path '*intermediate*' | head -1)
            cp "$ev" {output.events}
        fi

        # narrowPeak/bed (required)
        if [ -f "${{base}}.bed" ]; then
            cp "${{base}}.bed" {output.bed}
        else
            bd=$(find {params.outdir} -name '*.bed' -not -path '*intermediate*' | head -1)
            cp "$bd" {output.bed}
        fi

        # replicate-consistency table (may be absent if single replicate; touch if so)
        rc=$(find {params.outdir} -name '*replicationcodes*' -not -path '*intermediate*' | head -1 || true)
        if [ -n "${{rc:-}}" ] && [ -f "$rc" ]; then
            cp "$rc" {output.repcodes}
        else
            echo "# no replicationcodes table produced" > {output.repcodes}
        fi
        """


rule peaks:
    """Aggregate target: call peaks for every strain."""
    input:
        expand(os.path.join(PEAKS, "{strain}", "{strain}.events"), strain=STRAINS),
        expand(os.path.join(PEAKS, "{strain}", "{strain}.bed"), strain=STRAINS),
