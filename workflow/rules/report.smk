# =============================================================================
# report.smk — aggregate all QC into one MultiQC report + a decision scorecard
# =============================================================================

REPORT = os.path.join(RESULTS, "02_report")


def multiqc_inputs():
    ins = []
    for s in SAMPLE_IDS:
        ins += [
            os.path.join(QC, "samtools", f"{s}.flagstat.txt"),
            os.path.join(QC, "samtools", f"{s}.idxstats.txt"),
            os.path.join(QC, "samtools", f"{s}.stats.txt"),
            os.path.join(QC, "fastqc", f"{s}_R1_fastqc.zip"),
            os.path.join(QC, "fastqc", f"{s}_R2_fastqc.zip"),
        ]
    return ins


rule multiqc:
    input:
        multiqc_inputs(),
    output:
        html=os.path.join(REPORT, "multiqc_report.html"),
        data=directory(os.path.join(REPORT, "multiqc_report_data")),
    params:
        indir=QC,
        outdir=REPORT,
    log:
        os.path.join(LOG, "multiqc.log"),
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        multiqc {params.indir} -o {params.outdir} -n multiqc_report \
            --force > {log} 2>&1
        """


rule qc_scorecard:
    """Collate ChIPexoQual + samtools + header checks into one pass/flag table."""
    input:
        header=expand(os.path.join(QC, "header_check", "{s}.header_check.tsv"),
                      s=SAMPLE_IDS),
        flagstat=expand(os.path.join(QC, "samtools", "{s}.flagstat.txt"),
                        s=SAMPLE_IDS),
        chipexoqual=expand(
            os.path.join(QC, "chipexoqual", "{s}", "{s}.metrics.tsv"),
            s=SAMPLE_IDS),
    output:
        tsv=os.path.join(REPORT, "qc_scorecard.tsv"),
        md=os.path.join(REPORT, "qc_scorecard.md"),
    params:
        samples=SAMPLE_IDS,
        meta={s: SAMPLES[s] for s in SAMPLE_IDS},
        min_uniquely_mapped=config["qc"]["min_uniquely_mapped"],
        fsr_low=config["qc"]["fsr_balanced_low"],
        fsr_high=config["qc"]["fsr_balanced_high"],
    log:
        os.path.join(LOG, "scorecard.log"),
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/summarize_qc.py"
