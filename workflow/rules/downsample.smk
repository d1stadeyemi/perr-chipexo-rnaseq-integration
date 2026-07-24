# =============================================================================
# downsample.smk — depth-matched (symmetric) ChExMix control
# =============================================================================
# WHY THIS EXISTS
#   Primary calling gave A909 = 80 events, CJB111 = 28. But CJB111 was prepared
#   with half the standard bead/antibody slurry and sequenced shallower:
#     per-rep IP reads   A909 7.6-9.1M   vs  CJB111 3.4-4.9M
#     SignalFraction     A909 ~0.185     vs  CJB111 ~0.102
#   So the event gap confounds DEPTH with IP EFFICIENCY. This module removes the
#   depth term by subsampling EVERY replicate (both strains) to a single common
#   read count, and the controls likewise, then re-calling with identical
#   parameters.
#
# WHAT IT DOES AND DOES NOT SHOW
#   It isolates depth. It CANNOT equalise SignalFraction — that IS the IP
#   efficiency difference under test. So:
#     * if A909 collapses toward CJB111's count -> depth explains the gap
#     * if A909 stays high -> the residual gap is IP efficiency (and/or biology),
#       bounded but not eliminated. State this explicitly in the manuscript.
#
# Writes to a PARALLEL tree (results/04_downsample/) so primary peaks in
# results/03_peaks/ are never touched.
# =============================================================================

DS = os.path.join(RESULTS, "04_downsample")


rule ds_count_reads:
    """Primary-aligned, non-supplementary read count per sample (the ChExMix input)."""
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        cnt=os.path.join(DS, "counts", "{sample}.count"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {threads} -c -F 0x904 {input.bam} > {output.cnt}"


def all_ip_counts():
    return [os.path.join(DS, "counts", s + ".count") for s in IP_IDS]

def all_ctrl_counts():
    return [os.path.join(DS, "counts", s + ".count") for s in CTRL_IDS]


rule ds_targets:
    """Common target depths = min across all IP reps, and min across all controls."""
    input:
        ip=all_ip_counts(),
        ctrl=all_ctrl_counts(),
    output:
        tsv=os.path.join(DS, "targets.tsv"),
    run:
        def read(p):
            with open(p) as fh:
                return int(fh.read().strip())
        ip_counts = {os.path.basename(p)[:-6]: read(p) for p in input.ip}
        ct_counts = {os.path.basename(p)[:-6]: read(p) for p in input.ctrl}
        ip_target = min(ip_counts.values())
        ct_target = min(ct_counts.values())
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(output.tsv, "w") as fh:
            fh.write("kind\ttarget_reads\n")
            fh.write(f"IP\t{ip_target}\n")
            fh.write(f"control\t{ct_target}\n")
            fh.write("#\n#sample\tkind\tobserved\tfraction\n")
            for s, n in sorted(ip_counts.items()):
                fh.write(f"#{s}\tIP\t{n}\t{ip_target/n:.6f}\n")
            for s, n in sorted(ct_counts.items()):
                fh.write(f"#{s}\tcontrol\t{n}\t{ct_target/n:.6f}\n")


rule ds_subsample:
    """Subsample one BAM to the common target for its kind (IP or control)."""
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
        cnt=os.path.join(DS, "counts", "{sample}.count"),
        targets=os.path.join(DS, "targets.tsv"),
    output:
        bam=os.path.join(DS, "bam", "{sample}.ds.bam"),
        bai=os.path.join(DS, "bam", "{sample}.ds.bam.bai"),
    params:
        kind=lambda wc: SAMPLES[wc.sample]["kind"],
        seed=config["downsample"]["seed"],
    log:
        os.path.join(RESULTS, "logs", "ds_subsample_{sample}.log"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        n=$(cat {input.cnt})
        t=$(awk -v k={params.kind} '$1==k {{print $2}}' {input.targets})
        frac=$(awk -v t="$t" -v n="$n" 'BEGIN{{f=t/n; print (f>=1)?1:f}}')
        echo "sample={wildcards.sample} kind={params.kind} n=$n target=$t frac=$frac" > {log}
        if [ "$frac" = "1" ]; then
            cp {input.bam} {output.bam}
        else
            samtools view -@ {threads} -b \
                --subsample "$frac" --subsample-seed {params.seed} \
                {input.bam} > {output.bam}
        fi
        samtools index -@ {threads} {output.bam}
        samtools view -@ {threads} -c -F 0x904 {output.bam} >> {log}
        """


# --- ChExMix on the depth-matched BAMs (mirrors rule chexmix exactly) ---------
def ds_bam(sample):
    return os.path.join(DS, "bam", sample + ".ds.bam")

def ds_bai(sample):
    return os.path.join(DS, "bam", sample + ".ds.bam.bai")

def ds_signal_flags(strain):
    return " ".join(
        f"--expt{strain}-rep{i} {ds_bam(s)}"
        for i, s in enumerate(strain_ids(strain, "IP"), start=1)
    )

def ds_control_flags(strain):
    ctrl = strain_ids(strain, "control")
    if not ctrl:
        return ""
    cb = ds_bam(ctrl[0])
    n_ip = len(strain_ids(strain, "IP"))
    return " ".join(f"--ctrl{strain}-rep{i} {cb}" for i in range(1, n_ip + 1))

def ds_inputs(strain):
    ids = strain_ids(strain)
    return [ds_bam(s) for s in ids] + [ds_bai(s) for s in ids]


rule ds_chexmix:
    """Re-call peaks on depth-matched data.

    CRITICAL: every ChExMix parameter here is read from the SAME config keys as
    rule `chexmix` (mappability, q, minfold, format, exclude). If the control
    were called with different parameters than the primary run, the comparison
    would be meaningless. Do not hard-code values in this rule."""
    input:
        bams=lambda wc: ds_inputs(wc.strain),
        geninfo=lambda wc: ref_geninfo(strain_assembly(wc.strain)),
    output:
        events=os.path.join(DS, "peaks", "{strain}_ds", "{strain}_ds.events"),
    params:
        outdir=os.path.join(DS, "peaks", "{strain}_ds"),
        outname=lambda wc: f"{wc.strain}_ds",
        signal=lambda wc: ds_signal_flags(wc.strain),
        control=lambda wc: ds_control_flags(wc.strain),
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
        os.path.join(RESULTS, "logs", "ds_chexmix_{strain}.log"),
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

        # ChExMix writes into OUTDIR/NAME/ with files prefixed NAME_NAME.
        base="{params.outdir}/{params.outname}/{params.outname}_{params.outname}"
        if [ -f "${{base}}.events" ]; then
            cp "${{base}}.events" {output.events}
        else
            ev=$(find {params.outdir} -name '*.events' -not -path '*intermediate*' | head -1)
            cp "$ev" {output.events}
        fi
        """


rule ds_summary:
    """Side-by-side: primary vs depth-matched event counts per strain."""
    input:
        primary=expand(os.path.join(PEAKS, "{strain}", "{strain}.events"), strain=STRAINS),
        ds=expand(os.path.join(DS, "peaks", "{strain}_ds", "{strain}_ds.events"), strain=STRAINS),
        targets=os.path.join(DS, "targets.tsv"),
    output:
        tsv=os.path.join(DS, "downsample_summary.tsv"),
    run:
        def n_events(p):
            with open(p) as fh:
                return sum(1 for ln in fh if not ln.startswith("#"))
        rows = []
        for strain in STRAINS:
            p = os.path.join(PEAKS, strain, f"{strain}.events")
            d = os.path.join(DS, "peaks", f"{strain}_ds", f"{strain}_ds.events")
            rows.append((strain, n_events(p), n_events(d)))
        with open(output.tsv, "w") as fh:
            fh.write("strain\tevents_primary\tevents_depth_matched\n")
            for r in rows:
                fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")


rule downsample:
    """Aggregate target."""
    input:
        os.path.join(DS, "downsample_summary.tsv"),
