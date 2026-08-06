# =============================================================================
# concordance.smk — cross-caller agreement between ChExMix and GEM
# =============================================================================
# WHY THIS IS A SEPARATE, MINIMAL MODULE
#   GEM is NOT wired in as a peak caller, deliberately. Making it run on this
#   dataset required hand-calibrating two parameters (--icr and --a) because
#   three of its estimators degrade at this read depth: on a 2.1 Mb genome with
#   76% of positions covered, GEM's control-scaling estimator collapsed to 0.01,
#   its Poisson expectation reached 8,352 reads per 301 bp window, and its
#   auto-estimated sparse prior admitted 2 events from 56 M reads.
#
#   A Snakemake rule hard-coding `--a 50 --icr 5.88` would look principled and
#   would not be: those values are specific to these BAMs, and anyone running
#   the workflow on other data would get nonsense with no signal that
#   recalibration is needed. The GEM run therefore lives in a documented script
#   (workflow/scripts/gem_concordance.sh) whose header states what was tuned and
#   why.
#
#   What DOES belong in the pipeline is the NUMBER. This module recomputes the
#   overlap from files already on disk, so the concordance figure quoted in the
#   manuscript regenerates from the workflow even though the GEM run does not.
#
# WHAT IT MEASURES
#   For each tolerance, how many GEM events fall within that distance of a
#   ChExMix peak, and vice versa. Reported at several tolerances because the
#   right one is not obvious: GEM resolves the perR locus into ~5 components
#   where ChExMix reports 2, so a strict tolerance understates agreement at the
#   locus level while a loose one overstates it at the event level.
# =============================================================================

CONC = os.path.join(RESULTS, "12_concordance")


rule concordance_table:
    """Overlap between ChExMix and GEM event sets, at several tolerances."""
    input:
        chexmix=os.path.join(PEAKS_CHEXMIX, "{strain}", "{strain}.events"),
        gem=lambda wc: config["concordance"]["gem_events"][wc.strain],
    output:
        tsv=os.path.join(CONC, "{strain}.concordance.tsv"),
        detail=os.path.join(CONC, "{strain}.event_pairs.tsv"),
    params:
        strain=lambda wc: wc.strain,
        tolerances=config["concordance"]["tolerances"],
        flagship=config["concordance"].get("flagship", {}),
    log:
        os.path.join(LOG, "concordance_{strain}.log"),
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/concordance.py"


rule concordance:
    """Aggregate target."""
    input:
        expand(os.path.join(CONC, "{strain}.concordance.tsv"),
               strain=[s for s in STRAINS
                       if s in config.get("concordance", {}).get("gem_events", {})]),
