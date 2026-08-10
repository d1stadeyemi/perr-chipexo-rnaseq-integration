# =============================================================================
# setup.smk — unpack the sequencing-core delivery
# =============================================================================
# The core's zip triggers unzip's "possible zip bomb" heuristic (overlapping
# components). We extract with Python's zipfile instead, which is not fooled by
# that heuristic, so the workflow works out-of-the-box with no env variables.
# =============================================================================


# Extraction applies only to the archive entry point. In fastq mode the BAMs are
# produced by align.smk instead, and defining this rule unconditionally would
# leave two rules capable of producing results/00_raw/*.bam. Snakemake resolves
# that ambiguity silently in favour of whichever it sees first, which meant the
# alignment rules never ran.
if INPUT_MODE == "zip":

    rule stage_delivery:
        """Extract the delivery archive into results/00_raw/ (structure preserved)."""
        input:
            zip=INPUT_ZIP,
        output:
            # sentinel + every raw file the rest of the DAG consumes
            flag=touch(os.path.join(RAW, ".extracted")),
            files=all_raw_targets(),
        log:
            os.path.join(RESULTS, "logs", "stage_delivery.log"),
        params:
            outdir=RAW,
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/stage_delivery.py"
