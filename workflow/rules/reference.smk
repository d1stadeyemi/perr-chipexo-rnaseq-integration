# =============================================================================
# reference.smk — acquire and prepare reference genomes for peak calling
# =============================================================================
# INTERNET REQUIRED for `get_reference` (NCBI download). On HPCs where compute
# nodes are offline, run this stage on a login/transfer node first:
#     snakemake --sdm conda -c1 get_reference_all
# Downstream peak calling then runs offline on a compute/interactive node.
#
# NOTE: Snakemake output paths must be STATIC templates (no functions), so we
# write them with the {asm} wildcard directly; they mirror the ref_*() helpers
# in common.smk (REF/{asm}/{asm}.fna etc.).
# =============================================================================


rule get_reference:
    """Download genome FASTA + GFF3 for one accession via the NCBI datasets CLI."""
    output:
        fna=os.path.join(REF, "{asm}", "{asm}.fna"),
        gff=os.path.join(REF, "{asm}", "{asm}.gff"),
    params:
        asm=lambda wc: wc.asm,
        dldir=os.path.join(REF, "{asm}", "_download"),
        zipf=os.path.join(REF, "{asm}", "_dl.zip"),
    log:
        os.path.join(RESULTS, "logs", "get_reference_{asm}.log"),
    conda:
        "../envs/reference.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dldir}
        datasets download genome accession {params.asm} \
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
        fna=os.path.join(REF, "{asm}", "{asm}.fna"),
    output:
        fai=os.path.join(REF, "{asm}", "{asm}.fna.fai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.fna}"


rule reference_geninfo:
    """ChExMix --geninfo file: chrName<TAB>chrLength (derived from the .fai)."""
    input:
        fai=os.path.join(REF, "{asm}", "{asm}.fna.fai"),
    output:
        geninfo=os.path.join(REF, "{asm}", "{asm}.geninfo"),
    shell:
        "cut -f1,2 {input.fai} > {output.geninfo}"


rule reference_seqdir:
    """Per-chromosome FASTA directory for ChExMix --seq (used in motif stage)."""
    input:
        fna=os.path.join(REF, "{asm}", "{asm}.fna"),
    output:
        seqdir=directory(os.path.join(REF, "{asm}", "seq")),
    shell:
        r"""
        mkdir -p {output.seqdir}
        awk -v d={output.seqdir} \
            '/^>/{{name=substr($1,2); f=d"/"name".fa"}} {{print > f}}' {input.fna}
        """


rule get_reference_all:
    """Convenience target: fetch + prep every accession (run on a login node)."""
    input:
        [ref_geninfo(a) for a in all_assemblies()],
        [ref_seqdir(a)  for a in all_assemblies()],
        [ref_gff(a)     for a in all_assemblies()],
