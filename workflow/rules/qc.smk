# =============================================================================
# qc.smk — quality control of the core's pre-aligned ChIP-exo data
# =============================================================================
# We do NOT re-align here: the delivery is already mapped + deduplicated by the
# core. This stage independently certifies that the data is sound before any
# peak calling. Tools:
#   samtools     -> provenance (genome build, mapping, pairing) + MultiQC stats
#   FastQC       -> raw read quality / adapter content (short inserts expected)
#   ChIPexoQual  -> ChIP-exo-specific QC (ARC/URC enrichment, FSR strandedness,
#                   library complexity via beta1/beta2)   [Welch et al. 2017]
#   deepTools    -> IP-vs-control enrichment (fingerprint), replicate
#                   correlation, PE fragment-size distribution
# =============================================================================

QC = os.path.join(RESULTS, "01_qc")
LOG = os.path.join(RESULTS, "logs")


# ---------------------------------------------------------------------------
# BAM indexing + header/provenance verification
# ---------------------------------------------------------------------------
rule samtools_index:
    input:
        bam=lambda wc: raw_bam(wc.sample),
        flag=os.path.join(RAW, ".extracted"),
    output:
        bai=os.path.join(RAW, "844_LBurcham_UTenn", "{sample}" + BAM_SUF + ".bai"),
    log:
        os.path.join(LOG, "index_{sample}.log"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        # index the ORIGINAL bam in place next to it, then expose via output path
        "samtools index -@ {threads} {input.bam} {output.bai} 2> {log}"


rule check_bam_header:
    """Fail loudly if the @SQ genome build != the accession declared in config."""
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        report=os.path.join(QC, "header_check", "{sample}.header_check.tsv"),
    params:
        accession=lambda wc: SAMPLES[wc.sample]["accession"],
        chrom_len=lambda wc: SAMPLES[wc.sample]["chrom_len"],
        sample_id=lambda wc: wc.sample,
    log:
        os.path.join(LOG, "header_{sample}.log"),
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/check_bam_header.py"


rule samtools_flagstat:
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        os.path.join(QC, "samtools", "{sample}.flagstat.txt"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"


rule samtools_idxstats:
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        os.path.join(QC, "samtools", "{sample}.idxstats.txt"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools idxstats {input.bam} > {output}"


rule samtools_stats:
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        os.path.join(QC, "samtools", "{sample}.stats.txt"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats -@ {threads} {input.bam} > {output}"


# ---------------------------------------------------------------------------
# FastQC on the raw reads (both mates)
# ---------------------------------------------------------------------------
rule fastqc:
    input:
        fq=lambda wc: (raw_r1(wc.sample) if wc.mate == "R1" else raw_r2(wc.sample)),
        flag=os.path.join(RAW, ".extracted"),
    output:
        html=os.path.join(QC, "fastqc", "{sample}_{mate}_fastqc.html"),
        zip=os.path.join(QC, "fastqc", "{sample}_{mate}_fastqc.zip"),
    params:
        outdir=os.path.join(QC, "fastqc"),
        kmers=config["qc"]["fastqc_kmers"],
    log:
        os.path.join(LOG, "fastqc_{sample}_{mate}.log"),
    threads: 2
    conda:
        "../envs/fastqc.yaml"
    shell:
        # FastQC names outputs after the input basename; rename to our scheme.
        r"""
        set -euo pipefail
        tmp=$(mktemp -d)
        fastqc -k {params.kmers} -t {threads} -o "$tmp" {input.fq} > {log} 2>&1
        html=$(ls "$tmp"/*_fastqc.html)
        zipf=$(ls "$tmp"/*_fastqc.zip)
        mv "$html" {output.html}
        mv "$zipf" {output.zip}
        rm -rf "$tmp"
        """


# ---------------------------------------------------------------------------
# ChIPexoQual — the ChIP-exo-specific QC (per sample)
# ---------------------------------------------------------------------------

rule subsample_bam_for_qc:
    """Subsample to a fixed read count for ChIPexoQual island formation."""
    input:
        bam=lambda wc: raw_bam(wc.sample),
        bai=lambda wc: raw_bai(wc.sample),
    output:
        bam=os.path.join(QC, "subsampled", "{sample}.sub.bam"),
        bai=os.path.join(QC, "subsampled", "{sample}.sub.bam.bai"),
    params:
        target=config["qc"]["chipexoqual_target_reads"],
        seed=config["qc"]["chipexoqual_subsample_seed"],
    log:
        os.path.join(LOG, "subsample_{sample}.log"),
    threads: 4
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        set -euo pipefail
        total=$(samtools view -@ {threads} -c -F 0x904 {input.bam})
        frac=$(awk -v t={params.target} -v n="$total" 'BEGIN{{f=t/n; print (f>=1)?1:f}}')
        echo "total=$total target={params.target} frac=$frac" > {log}
        if [ "$frac" = "1" ]; then
            cp {input.bam} {output.bam}
        else
            samtools view -@ {threads} -b \
                --subsample "$frac" --subsample-seed {params.seed} \
                {input.bam} > {output.bam}
        fi
        samtools index -@ {threads} {output.bam}
        """

rule chipexoqual:
    input:
        bam=os.path.join(QC, "subsampled", "{sample}.sub.bam"),
        bai=os.path.join(QC, "subsampled", "{sample}.sub.bam.bai"),
    output:
        arc_urc=os.path.join(QC, "chipexoqual", "{sample}", "{sample}.ARCvURC.pdf"),
        fsr=os.path.join(QC, "chipexoqual", "{sample}", "{sample}.FSRdist.pdf"),
        regioncomp=os.path.join(QC, "chipexoqual", "{sample}", "{sample}.regionComp.pdf"),
        paramdist=os.path.join(QC, "chipexoqual", "{sample}", "{sample}.paramDist.pdf"),
        stats=os.path.join(QC, "chipexoqual", "{sample}", "{sample}.metrics.tsv"),
    params:
        sample_id=lambda wc: wc.sample,
        kind=lambda wc: SAMPLES[wc.sample]["kind"],
        strain=lambda wc: SAMPLES[wc.sample]["strain"],
        mc_cores=config["qc"]["chipexoqual_mc_cores"],
    log:
        os.path.join(LOG, "chipexoqual_{sample}.log"),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        "../envs/chipexoqual.yaml"
    script:
        "../scripts/run_chipexoqual.R"


# ---------------------------------------------------------------------------
# deepTools — enrichment fingerprint, replicate correlation, fragment size
# (grouped per strain: the 3 IP reps + the 1 no-tag control)
# ---------------------------------------------------------------------------
def strain_bams(strain):
    return [raw_bam(s) for s in strain_ids(strain)]

def strain_bais(strain):
    return [raw_bai(s) for s in strain_ids(strain)]

def strain_labels(strain):
    labs = []
    for s in strain_ids(strain):
        r = SAMPLES[s]
        tag = "IP" if r["kind"] == "IP" else "ctrl"
        labs.append(f"{r['strain']}_{tag}_{r['replicate']}")
    return labs


rule deeptools_fingerprint:
    """IP-vs-control enrichment: a good IP shows a strongly skewed fingerprint."""
    input:
        bams=lambda wc: strain_bams(wc.strain),
        bais=lambda wc: strain_bais(wc.strain),
    output:
        plot=os.path.join(QC, "deeptools", "{strain}.fingerprint.pdf"),
        counts=os.path.join(QC, "deeptools", "{strain}.fingerprint.counts.tsv"),
        metrics=os.path.join(QC, "deeptools", "{strain}.fingerprint.metrics.tsv"),
    params:
        labels=lambda wc: " ".join(strain_labels(wc.strain)),
        extra=config["qc"]["deeptools_extra"],
    log:
        os.path.join(LOG, "fingerprint_{strain}.log"),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=180,
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        plotFingerprint -b {input.bams} --labels {params.labels} \
            --plotFile {output.plot} \
            --outRawCounts {output.counts} \
            --outQualityMetrics {output.metrics} \
            --numberOfProcessors {threads} {params.extra} > {log} 2>&1
        """


rule deeptools_multibamsummary:
    input:
        bams=lambda wc: strain_bams(wc.strain),
        bais=lambda wc: strain_bais(wc.strain),
    output:
        npz=os.path.join(QC, "deeptools", "{strain}.multibam.npz"),
    params:
        labels=lambda wc: " ".join(strain_labels(wc.strain)),
        binsize=config["qc"]["deeptools_binsize"],
        extra=config["qc"]["deeptools_extra"],
    log:
        os.path.join(LOG, "multibam_{strain}.log"),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=180,
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        multiBamSummary bins --bamfiles {input.bams} \
            --labels {params.labels} --binSize {params.binsize} \
            --numberOfProcessors {threads} {params.extra} \
            -o {output.npz} > {log} 2>&1
        """


rule deeptools_correlation:
    """Spearman correlation heatmap: IP reps should cluster, control apart."""
    input:
        npz=os.path.join(QC, "deeptools", "{strain}.multibam.npz"),
    output:
        plot=os.path.join(QC, "deeptools", "{strain}.spearman.pdf"),
        matrix=os.path.join(QC, "deeptools", "{strain}.spearman.tsv"),
    log:
        os.path.join(LOG, "corr_{strain}.log"),
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        plotCorrelation -in {input.npz} \
            --corMethod spearman --whatToPlot heatmap \
            --skipZeros --plotNumbers \
            --plotFile {output.plot} \
            --outFileCorMatrix {output.matrix} > {log} 2>&1
        """


rule deeptools_fragsize:
    """Paired-end insert-size distribution (report noted ~67-85 bp average)."""
    input:
        bams=lambda wc: strain_bams(wc.strain),
        bais=lambda wc: strain_bais(wc.strain),
    output:
        plot=os.path.join(QC, "deeptools", "{strain}.fragsize.pdf"),
        table=os.path.join(QC, "deeptools", "{strain}.fragsize.tsv"),
    params:
        labels=lambda wc: " ".join(strain_labels(wc.strain)),
    log:
        os.path.join(LOG, "fragsize_{strain}.log"),
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    shell:
        r"""
        bamPEFragmentSize -b {input.bams} --samplesLabel {params.labels} \
            --numberOfProcessors {threads} \
            --maxFragmentLength 500 \
            --histogram {output.plot} \
            --outRawFragmentLengths {output.table} > {log} 2>&1
        """
