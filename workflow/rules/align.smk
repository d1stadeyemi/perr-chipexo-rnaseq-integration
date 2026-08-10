# =============================================================================
# align.smk -- optional entry point: align from FASTQ instead of using the
#              sequencing core's delivered BAMs
# =============================================================================
# WHY THIS MODULE EXISTS
#   The main workflow starts from the aligned BAMs supplied by the sequencing
#   core. That is the fastest path and it is what produced the published
#   results, but it is not what will be publicly archived: SRA holds FASTQ.
#   Anyone reproducing this work from the public record therefore needs to
#   align first, and needs to do it the same way the core did.
#
# WHAT THE CORE DID
#   Recovered verbatim from the @PG lines of the delivered BAMs:
#
#       bwa 0.7.17-r1188  ->  Picard MarkDuplicates 2.18.2  ->  samtools 1.21
#
#   Two details from flagstat that this module reproduces deliberately:
#     * duplicates are REMOVED, not marked (the delivered BAMs report 0
#       duplicates, so REMOVE_DUPLICATES=true, not the Picard default)
#     * no unmapped reads are retained
#   Getting either wrong changes every downstream count.
#
# HOW IT PLUGS IN
#   The rest of the workflow addresses samples by a CEGR-convention basename,
#   for example
#       43893_6xHis_ab9108_A9091_-_-_-_BX_FilteredBAM.bam
#   and parses strain, IP-vs-control and replicate number out of that name.
#   This module therefore writes its output to exactly those paths. Downstream
#   rules cannot tell the difference between a BAM produced here and one
#   unpacked from the core's archive.
#
#   FASTQ files carry no such metadata (an SRA accession is just SRR#######),
#   so the mapping from FASTQ to CEGR name is supplied by a small sample sheet.
#   See config/samples_fastq.tsv.example.
#
# ENABLING IT
#   config.yaml:
#       input_mode: "fastq"          # default is "zip"
#       alignment:
#         sample_sheet: "config/samples_fastq.tsv"
#         fastq_dir: "data/fastq"
#
#   In "zip" mode this module defines rules but nothing requests them, so it is
#   inert. The two entry points are mutually exclusive by construction: both
#   write results/00_raw/{sample}.bam, and Snakemake will refuse to build a DAG
#   in which two rules can produce the same file. Set input_mode to select one.
#
# NOT VALIDATED AGAINST THE PUBLISHED RUN
#   The published results were produced from the core's BAMs. This module
#   reproduces the core's documented procedure but has not itself been used to
#   regenerate them. Anyone using it should expect small differences from
#   aligner non-determinism (bwa mem is deterministic for a fixed thread count,
#   but thread count affects read ordering, which affects duplicate selection).
#   Verify with the header_check and flagstat rules before drawing conclusions.
# =============================================================================

ALIGN = os.path.join(RESULTS, "00_align")


# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------
# Sample discovery lives in common.smk, which is included first and handles both
# entry points. This module simply consumes the result, so there is exactly one
# definition of what a sample is.
FASTQ_SAMPLES = SAMPLES if config.get("input_mode") == "fastq" else {}


def _fq_strain_assembly(sample):
    return config["strains"][FASTQ_SAMPLES[sample]["strain"]]["assembly"]


# ---------------------------------------------------------------------------
# 1. BWA index
# ---------------------------------------------------------------------------
rule bwa_index:
    """Index the reference once per assembly.

    The .bwt file is the target; bwa writes five files with the FASTA as their
    common prefix and this rule declares all of them so a partial index cannot
    be mistaken for a complete one."""
    input:
        fna=os.path.join(REF, "{asm}", "{asm}.fna"),
    output:
        multiext(os.path.join(REF, "{asm}", "{asm}.fna"),
                 ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        os.path.join(LOG, "bwa_index_{asm}.log"),
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/align.yaml"
    shell:
        "bwa index {input.fna} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 2. Align, sort
# ---------------------------------------------------------------------------
rule bwa_mem:
    """Paired-end alignment, coordinate-sorted.

    -R writes a read group. Picard MarkDuplicates requires one and fails
    without it, which is a common and confusing first error when adapting this
    module. SM is set to the CEGR sample id so the read group is traceable.

    Sorting is piped rather than written to an intermediate SAM: on a 2 Mb
    genome the alignment is small, but the pipe avoids a large temporary file
    for no cost."""
    input:
        r1=lambda wc: raw_r1(wc.sample),
        r2=lambda wc: raw_r2(wc.sample),
        idx=lambda wc: multiext(
            os.path.join(REF, _fq_strain_assembly(wc.sample),
                         f"{_fq_strain_assembly(wc.sample)}.fna"),
            ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam=temp(os.path.join(ALIGN, "{sample}.sorted.bam")),
    params:
        fna=lambda wc: os.path.join(
            REF, _fq_strain_assembly(wc.sample),
            f"{_fq_strain_assembly(wc.sample)}.fna"),
        rg=lambda wc: (
            r"@RG\tID:" + wc.sample + r"\tSM:" + wc.sample +
            r"\tPL:ILLUMINA\tLB:" + wc.sample),
    log:
        os.path.join(LOG, "bwa_mem_{sample}.log"),
    threads: config.get("alignment", {}).get("bwa_threads", 8)
    resources:
        mem_mb=16000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        bwa mem -t {threads} -R '{params.rg}' \
            {params.fna} {input.r1} {input.r2} 2> {log} \
          | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        """


# ---------------------------------------------------------------------------
# 3. Remove duplicates, filter, index
# ---------------------------------------------------------------------------
rule mark_duplicates:
    """Picard MarkDuplicates with REMOVE_DUPLICATES=true, then drop unmapped
    reads, producing a BAM equivalent to the core's FilteredBAM.

    REMOVE_DUPLICATES=true is deliberate and is NOT the Picard default. The
    delivered BAMs report zero duplicates in flagstat, which means the core
    removed them rather than flagging them. Leaving the default would carry
    duplicate reads into every downstream count.

    -F 0x904 excludes unmapped (0x4), secondary (0x100) and supplementary
    (0x800) alignments. This is the same filter the QC module uses when
    reporting primary mapped read counts, so the two agree by construction."""
    input:
        bam=os.path.join(ALIGN, "{sample}.sorted.bam"),
    output:
        bam=os.path.join(RAW, "{sample}" + config["naming"]["bam_suffix"]),
        bai=os.path.join(RAW, "{sample}" + config["naming"]["bam_suffix"] + ".bai"),
        metrics=os.path.join(ALIGN, "{sample}.dupmetrics.txt"),
    params:
        mem=config.get("alignment", {}).get("picard_mem_gb", 8),
        tmp=os.path.join(ALIGN, "tmp_{sample}"),
    log:
        os.path.join(LOG, "markdup_{sample}.log"),
    threads: 4
    resources:
        mem_mb=lambda wc: (config.get("alignment", {}).get("picard_mem_gb", 8) + 4) * 1000,
        runtime=120,
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.tmp} $(dirname {output.bam})
        picard -Xmx{params.mem}g MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={params.tmp}/dedup.bam \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR={params.tmp} > {log} 2>&1
        samtools view -@ {threads} -b -F 0x904 \
            {params.tmp}/dedup.bam > {output.bam} 2>> {log}
        samtools index -@ {threads} {output.bam}
        rm -rf {params.tmp}
        echo "primary mapped: $(samtools view -c -F 0x904 {output.bam})" >> {log}
        """


# ---------------------------------------------------------------------------
# 4. Provenance record
# ---------------------------------------------------------------------------
rule align_manifest:
    """Record what was aligned and with what, so a FASTQ-derived run can be
    told apart from one built on the core's BAMs."""
    input:
        bams=lambda wc: expand(
            os.path.join(RAW, "{s}" + config["naming"]["bam_suffix"]),
            s=sorted(FASTQ_SAMPLES)),
    output:
        tsv=os.path.join(ALIGN, "alignment_manifest.tsv"),
    params:
        samples=lambda wc: FASTQ_SAMPLES,
    log:
        os.path.join(LOG, "align_manifest.log"),
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        printf 'sample_id\tprimary_mapped\tbwa_version\tpicard_version\n' > {output.tsv}
        bwa_v=$(bwa 2>&1 | grep -m1 Version | awk '{{print $2}}')
        pic_v=$(picard MarkDuplicates --version 2>&1 | head -1)
        for b in {input.bams}; do
            s=$(basename "$b" {config[naming][bam_suffix]})
            n=$(samtools view -c -F 0x904 "$b")
            printf '%s\t%s\t%s\t%s\n' "$s" "$n" "$bwa_v" "$pic_v" >> {output.tsv}
        done
        cat {output.tsv} > {log}
        """


rule align:
    """Aggregate target: align every sample in the sheet and record provenance."""
    input:
        os.path.join(ALIGN, "alignment_manifest.tsv"),
