# =============================================================================
# orthology.smk — cross-serotype comparison of the PerR regulon
# =============================================================================
# A909 uses SAK_RS* locus tags and CJB111 uses ID870_RS*, so the two regulon
# tables cannot be compared directly. This module builds 1:1 orthologs by
# RECIPROCAL BEST HITS between the two proteomes and re-expresses both regulons
# on that shared coordinate system.
#
# CAVEAT carried into the output: CJB111 yielded ~4x fewer binding events, and
# it was prepared with half the standard bead/antibody slurry at lower depth.
# "A909-only direct" is therefore partly a SENSITIVITY difference and must be
# read alongside the depth-matched control (results/04_downsample).
# =============================================================================

ORTH = os.path.join(RESULTS, "08_orthology")

# the two strains, in a fixed order
STRAIN_A = STRAINS[0]
STRAIN_B = STRAINS[1] if len(STRAINS) > 1 else STRAINS[0]


rule get_proteome:
    """Protein FASTA per assembly (NCBI datasets). Needs internet: run on a
    login/transfer node, e.g. snakemake -c1 get_proteome_all"""
    output:
        faa=os.path.join(REF, "{asm}", "{asm}.faa"),
    params:
        dldir=os.path.join(REF, "{asm}", "_prot"),
        zipf=os.path.join(REF, "{asm}", "_prot.zip"),
    log:
        os.path.join(RESULTS, "logs", "get_proteome_{asm}.log"),
    conda:
        "../envs/reference.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dldir}
        datasets download genome accession {wildcards.asm} \
            --include protein --filename {params.zipf} > {log} 2>&1
        unzip -o {params.zipf} -d {params.dldir} >> {log} 2>&1
        faa=$(find {params.dldir} -name 'protein.faa' | head -1)
        cp "$faa" {output.faa}
        rm -rf {params.dldir} {params.zipf}
        """


rule get_proteome_all:
    input:
        [os.path.join(REF, a, f"{a}.faa") for a in all_assemblies()],


rule blastp_pair:
    """One direction of the reciprocal search (query -> subject)."""
    input:
        q=os.path.join(REF, "{qasm}", "{qasm}.faa"),
        s=os.path.join(REF, "{sasm}", "{sasm}.faa"),
    output:
        tsv=os.path.join(ORTH, "blast", "{qasm}__vs__{sasm}.tsv"),
    params:
        evalue=config["orthology"]["evalue"],
        dbdir=os.path.join(ORTH, "blast", "db_{sasm}"),
    log:
        os.path.join(RESULTS, "logs", "blastp_{qasm}__vs__{sasm}.log"),
    threads: 4
    conda:
        "../envs/blast.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dbdir}
        makeblastdb -in {input.s} -dbtype prot -out {params.dbdir}/db > {log} 2>&1
        blastp -query {input.q} -db {params.dbdir}/db \
            -evalue {params.evalue} -max_target_seqs 5 -num_threads {threads} \
            -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
            -out {output.tsv} >> {log} 2>&1
        """


rule ortholog_table:
    """Reciprocal best hits -> 1:1 locus-tag ortholog pairs."""
    input:
        fwd=os.path.join(ORTH, "blast",
                         f"{strain_assembly(STRAIN_A)}__vs__{strain_assembly(STRAIN_B)}.tsv"),
        rev=os.path.join(ORTH, "blast",
                         f"{strain_assembly(STRAIN_B)}__vs__{strain_assembly(STRAIN_A)}.tsv"),
        gff_a=os.path.join(REF, strain_assembly(STRAIN_A),
                           f"{strain_assembly(STRAIN_A)}.gff"),
        gff_b=os.path.join(REF, strain_assembly(STRAIN_B),
                           f"{strain_assembly(STRAIN_B)}.gff"),
    output:
        tsv=os.path.join(ORTH, "orthologs_rbh.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
        min_identity=config["orthology"]["min_identity"],
        min_coverage=config["orthology"]["min_coverage"],
    log:
        os.path.join(RESULTS, "logs", "ortholog_table.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/build_ortholog_table.py"


rule crossstrain_regulon:
    """Re-express both regulons on the ortholog coordinate system."""
    input:
        reg_a=os.path.join(INTEG, STRAIN_A, f"{STRAIN_A}.direct_regulon.tsv"),
        reg_b=os.path.join(INTEG, STRAIN_B, f"{STRAIN_B}.direct_regulon.tsv"),
        orthologs=os.path.join(ORTH, "orthologs_rbh.tsv"),
    output:
        tsv=os.path.join(ORTH, "crossstrain_regulon.tsv"),
        summary=os.path.join(ORTH, "crossstrain_summary.tsv"),
    params:
        strain_a=STRAIN_A,
        strain_b=STRAIN_B,
    log:
        os.path.join(RESULTS, "logs", "crossstrain_regulon.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/crossstrain_regulon.py"


rule orthology:
    """Aggregate target."""
    input:
        os.path.join(ORTH, "crossstrain_regulon.tsv"),
        os.path.join(ORTH, "crossstrain_summary.tsv"),
