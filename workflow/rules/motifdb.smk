# =============================================================================
# motifdb.smk — compare the discovered motif against curated motif databases
# =============================================================================
# WHY THIS EXISTS
#   Our earlier TOMTOM runs compared the discovered motif against a target
#   database containing exactly ONE motif (the GAS Per box). TOMTOM warned:
#     "Target database size too small (1) for accurate p-value computation.
#      Provide at least 50 motifs."
#   So those p/q-values are NOT publishable. Comparing against curated databases
#   with hundreds of motifs gives a properly estimated null.
#
#   DATABASE CHOICE (MEME Suite "Combined Prokaryotes", 558 motifs total):
#     prodoric_2021.9.meme  333 motifs, w 5-34 (avg 15.5) - largest, most recent
#     regtransbase.meme     141 motifs, w 14-38 (avg 19.8) - curated, older resource
#     collectf.meme          84 motifs, w 9-53 (avg 21.1) - highest curation standard
#   All three are used together: TOMTOM needs >=50 motifs for a valid null, and
#   558 estimates it far better than any single set. A larger null also means
#   stricter multiple-testing correction - that is correct, not a problem.
#
#   DELIBERATELY EXCLUDED: fan2020.meme (115 motifs) is Pseudomonas syringae
#   HT-SELEX only - phylogenetically distant from Streptococcus, so it would
#   inflate the null without adding relevant coverage.
#
#   The real question this asks: does the GBS PerR motif resemble any known
#   bacterial TF motif - a Fur-family site, or something else entirely?
#
# The download needs internet: run `get_motif_databases` on a login node.
# =============================================================================

MDB = os.path.join("resources", "motif_databases")


def motifdb_files():
    """Absolute paths of the prokaryotic databases selected in config."""
    return [os.path.join(MDB, "PROKARYOTE", f)
            for f in config["motifdb"]["databases"]]


rule get_motif_databases:
    """Download + extract the MEME Suite curated motif databases (INTERNET).

    NOTE: the archive contains every MEME database (JASPAR, CIS-BP, miRBase,
    ...) and is large; only the PROKARYOTE files are used here."""
    output:
        dbs=motifdb_files(),
    params:
        url=config["motifdb"]["url"],
        tgz=os.path.join(MDB, "_motif_databases.tgz"),
        outdir=MDB,
    log:
        os.path.join(RESULTS, "logs", "get_motif_databases.log"),
    conda:
        "../envs/reference.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        curl -L -o {params.tgz} {params.url} > {log} 2>&1
        # the tarball unpacks to motif_databases/... ; strip that leading dir
        tar -xzf {params.tgz} -C {params.outdir} --strip-components=1 >> {log} 2>&1
        rm -f {params.tgz}
        ls -la {params.outdir}/PROKARYOTE >> {log} 2>&1
        """


rule tomtom_vs_databases:
    """TOMTOM: discovered motifs vs all curated prokaryotic motifs."""
    input:
        # depend on the MEME-ChIP sentinel: combined.meme is written by that rule
        # but is not one of its declared outputs.
        memechip_done=os.path.join(MOTIF, "{strain}", ".memechip_{peakset}.done"),
        dbs=motifdb_files(),
    output:
        tsv=os.path.join(MOTIF, "{strain}", "tomtom_db_{peakset}", "tomtom.tsv"),
    params:
        outdir=os.path.join(MOTIF, "{strain}", "tomtom_db_{peakset}"),
        memechip_dir=os.path.join(MOTIF, "{strain}", "memechip_{peakset}"),
        thresh=config["motifdb"]["thresh"],
        dist=config["motifdb"]["dist"],
    log:
        os.path.join(RESULTS, "logs", "tomtom_db_{strain}_{peakset}.log"),
    conda:
        "../envs/meme.yaml"
    shell:
        r"""
        set -euo pipefail
        q="{params.memechip_dir}/combined.meme"
        if [ ! -f "$q" ]; then
            q=$(find {params.memechip_dir} -name '*.meme' -not -name 'GAS_*' | head -1)
        fi
        tomtom -oc {params.outdir} \
            -thresh {params.thresh} -dist {params.dist} \
            -eps -png \
            "$q" \
            {input.dbs} \
            > {log} 2>&1
        """


rule motifdb_summary:
    """Collate the best database matches per strain x peak set."""
    input:
        expand(os.path.join(MOTIF, "{strain}", "tomtom_db_{peakset}", "tomtom.tsv"),
               strain=STRAINS, peakset=PEAKSETS),
    output:
        tsv=os.path.join(MOTIF, "motif_database_matches.tsv"),
    params:
        strains=STRAINS,
        peaksets=PEAKSETS,
        motifdir=MOTIF,
    log:
        os.path.join(RESULTS, "logs", "motifdb_summary.log"),
    run:
        import os as _os
        rows = []
        for st in params.strains:
            for ps in params.peaksets:
                p = _os.path.join(params.motifdir, st, f"tomtom_db_{ps}", "tomtom.tsv")
                if not _os.path.exists(p):
                    continue
                with open(p) as fh:
                    hdr = fh.readline().rstrip("\n").split("\t")
                    idx = {c: i for i, c in enumerate(hdr)}
                    for ln in fh:
                        if ln.startswith("#") or not ln.strip():
                            continue
                        f = ln.rstrip("\n").split("\t")
                        if len(f) < len(hdr):
                            continue
                        rows.append([st, ps] + [f[idx.get(c, 0)] for c in
                                    ("Query_ID", "Target_ID", "Optimal_offset",
                                     "p-value", "E-value", "q-value",
                                     "Overlap", "Orientation")])
        with open(output.tsv, "w") as ob:
            ob.write("strain\tpeakset\tquery_motif\tdb_match\toffset\t"
                     "p_value\tE_value\tq_value\toverlap\torientation\n")
            for r in rows:
                ob.write("\t".join(r) + "\n")
        with open(log[0], "w") as lg:
            lg.write(f"database match rows: {len(rows)}\n")


rule motifdb_family_check:
    """List Fur-family motifs present in the databases.

    A 'no match' result is only meaningful if the databases actually CONTAIN
    Fur-family sites to match against. This rule makes that explicit."""
    input:
        dbs=motifdb_files(),
    output:
        tsv=os.path.join(MOTIF, "database_fur_family_motifs.tsv"),
    shell:
        r"""
        set -euo pipefail
        {{
          echo -e "database\tmotif_line"
          for f in {input.dbs}; do
            grep -i '^MOTIF' "$f" \
              | grep -iE 'perr|[^a-z]fur|zur|mntr|nur|irr|troR' \
              | sed "s|^|$(basename $f)\t|" || true
          done
        }} > {output.tsv}
        """


rule motifdb:
    """Aggregate target."""
    input:
        os.path.join(MOTIF, "motif_database_matches.tsv"),
        os.path.join(MOTIF, "database_fur_family_motifs.tsv"),
