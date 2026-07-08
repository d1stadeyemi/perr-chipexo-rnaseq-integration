# =============================================================================
# reference.smk — acquire and prepare reference genomes for peak calling
# =============================================================================
# INTERNET REQUIRED for `get_reference` (NCBI download). On HPCs where compute
# nodes are offline, run this stage on a login/transfer node first:
#     snakemake --sdm conda -c1 get_reference_all
# Downstream peak calling then runs offline on a compute/interactive node.
#
# NOTE: Snakemake output paths must be STATIC templates (no functions), so we
# write them with the {acc} wildcard directly; they mirror the ref_*() helpers
# in common.smk (REF/{acc}/{acc}.fna etc.).
# =============================================================================


rule get_reference:
    """Download genome FASTA + GFF3 for one accession via the NCBI datasets CLI."""
    output:
        fna=os.path.join(REF, "{acc}", "{acc}.fna"),
        gff=os.path.join(REF, "{acc}", "{acc}.gff"),
    params:
        acc=lambda wc: wc.acc,
        dldir=os.path.join(REF, "{acc}", "_download"),
        zipf=os.path.join(REF, "{acc}", "_dl.zip"),
    log:
        os.path.join(RESULTS, "logs", "get_reference_{acc}.log"),
    conda:
        "../envs/reference.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dldir}
        datasets download genome accession {params.acc} \
            --include genome,gff3 --filename {params.zipf} > {log} 2>&1
        unzip -o {params.zipf} -d {params.dldir} >> {log} 2>&1
        fna=$(find {params.dldir} -name '*.fna' | head -1)
        gff=$(find {params.dldir} -name 'genomic.gff' | head -1)
        cp "$fna" {output.fna}
        cp "$gff" {output.gff}
        rm -rf {params.dldir} {params.zipf}
        """


rule reference_faidx:
    input:
        fna=os.path.join(REF, "{acc}", "{acc}.fna"),
    output:
        fai=os.path.join(REF, "{acc}", "{acc}.fna.fai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.fna}"


rule reference_geninfo:
    """ChExMix --geninfo file: chrName<TAB>chrLength (derived from the .fai)."""
    input:
        fai=os.path.join(REF, "{acc}", "{acc}.fna.fai"),
    output:
        geninfo=os.path.join(REF, "{acc}", "{acc}.geninfo"),
    shell:
        "cut -f1,2 {input.fai} > {output.geninfo}"


rule reference_seqdir:
    """Per-chromosome FASTA directory for ChExMix --seq (used in motif stage)."""
    input:
        fna=os.path.join(REF, "{acc}", "{acc}.fna"),
    output:
        seqdir=directory(os.path.join(REF, "{acc}", "seq")),
    shell:
        r"""
        mkdir -p {output.seqdir}
        awk -v d={output.seqdir} \
            '/^>/{{name=substr($1,2); f=d"/"name".fa"}} {{print > f}}' {input.fna}
        """


rule get_reference_all:
    """Convenience target: fetch + prep every accession (run on a login node)."""
    input:
        [ref_geninfo(a) for a in all_accessions()],
        [ref_seqdir(a) for a in all_accessions()],
        [ref_gff(a) for a in all_accessions()],
