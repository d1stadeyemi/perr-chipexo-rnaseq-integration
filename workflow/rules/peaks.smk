# =============================================================================
# peaks.smk — ChExMix peak (binding-event) calling
# =============================================================================
# ChExMix (Yamada et al. 2019; Pugh/Mahony) is the primary caller: a ChIP-exo
# mixture model that is control-aware and reports, per binding event, tag counts
# in signal vs the no-tag control and a log2 fold. One run per strain uses all
# three PerR-His replicates against that strain's no-tag control.
#
# We run in event-calling mode (--geninfo + reads + control). Motif-defined
# subtypes (--seq/--back) are deferred to the dedicated motif stage, where the
# AT-rich GBS background is built explicitly.
#
# ChExMix reads BAM through its SAM loader (--format SAM). Replicates and their
# matched control are passed with the --exptCOND-REP / --ctrlCOND-REP flags.
# =============================================================================


def chexmix_signal_flags(strain):
    """--exptSTRAIN-repN <bam> for each IP replicate."""
    parts = []
    for i, s in enumerate(strain_ids(strain, "IP"), start=1):
        parts.append(f"--expt{strain}-rep{i} {raw_bam(s)}")
    return " ".join(parts)


def chexmix_control_flags(strain):
    """Map the single no-tag control to every signal replicate of the strain."""
    ctrl = strain_ids(strain, "control")
    if not ctrl:
        return ""
    ctrl_bam = raw_bam(ctrl[0])
    n_ip = len(strain_ids(strain, "IP"))
    return " ".join(f"--ctrl{strain}-rep{i} {ctrl_bam}" for i in range(1, n_ip + 1))


def chexmix_input_bams(strain):
    bams = [raw_bam(s) for s in strain_ids(strain)]
    bais = [raw_bai(s) for s in strain_ids(strain)]
    return bams + bais


rule chexmix:
    """Call PerR binding events for one strain with all IP reps + control."""
    input:
        bams=lambda wc: chexmix_input_bams(wc.strain),
        geninfo=lambda wc: ref_geninfo(strain_assembly(wc.strain)),
    output:
        events=os.path.join(PEAKS, "{strain}", "{strain}.events"),
        done=touch(os.path.join(PEAKS, "{strain}", ".chexmix.done")),
    params:
        outdir=os.path.join(PEAKS, "{strain}"),
        outname=lambda wc: wc.strain,
        signal=lambda wc: chexmix_signal_flags(wc.strain),
        control=lambda wc: chexmix_control_flags(wc.strain),
        fmt=config["peaks"]["format"],
        qval=config["peaks"]["chexmix_qvalue"],
        mem=config["peaks"]["chexmix_mem_gb"],
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
            --q {params.qval} \
            --out {params.outdir}/{params.outname} \
            {params.extra} > {log} 2>&1
        # ChExMix writes <out>_<...>.events; expose a stable path for the DAG
        ev=$(find {params.outdir} -name '*.events' | head -1)
        cp "$ev" {output.events}
        """

rule peaks:
    """Aggregate target: call peaks for every strain."""
    input:
        expand(os.path.join(PEAKS, "{strain}", "{strain}.events"), strain=STRAINS),
