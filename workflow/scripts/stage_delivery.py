"""Extract the sequencing-core delivery zip.

Uses Python's zipfile so the 'possible zip bomb (overlapped components)'
heuristic that `unzip` trips on is never an issue. Only files that downstream
rules actually need are guaranteed; everything else in the archive is extracted
too (report .pptx, sampleID.txt, bigWigs) for provenance.
"""

import os
import sys
import zipfile

# Snakemake injects `snakemake` into the namespace when run via `script:`.
zip_path = snakemake.input.zip                      # noqa: F821
outdir = snakemake.params.outdir                    # noqa: F821
logpath = snakemake.log[0]                           # noqa: F821

os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.dirname(logpath), exist_ok=True)

with open(logpath, "w") as log:
    log.write(f"Extracting {zip_path} -> {outdir}\n")
    with zipfile.ZipFile(zip_path) as zf:
        members = zf.namelist()
        log.write(f"{len(members)} entries in archive\n")
        for m in members:
            zf.extract(m, path=outdir)
            log.write(f"  extracted {m}\n")
    log.write("done\n")

# sanity: confirm the explicit downstream targets now exist
missing = [f for f in snakemake.output.files if not os.path.exists(f)]  # noqa: F821
if missing:
    sys.exit("ERROR: expected files missing after extraction:\n  "
             + "\n  ".join(missing))
