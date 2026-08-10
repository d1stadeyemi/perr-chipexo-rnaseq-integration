# =============================================================================
# common.smk — sample discovery, helpers, global wildcards
# =============================================================================
# Sample discovery happens at DAG-construction time by reading the ZIP central
# directory (zipfile.namelist()) WITHOUT extracting. This is what lets the user
# "just hand over the zip": we know every sample before a single byte is
# unpacked. Extraction itself is a normal rule (rule stage_delivery).
# =============================================================================

import os
import re
import zipfile
from pathlib import PurePosixPath


configfile: "config/config.yaml"


RESULTS = config["results_dir"]
RAW = os.path.join(RESULTS, "00_raw")
# Defined here rather than in qc.smk so that every module can use it regardless
# of include order.
LOG = os.path.join(RESULTS, "logs")
# Only required in zip mode. Read with .get() so a fastq-only user does not
# need to declare an archive they do not have.
INPUT_ZIP = config.get("input_zip")

# --- map each strain PREFIX -> (strain_name, accession, chrom_len) ------------
_STRAIN_BY_PREFIX = {}
for sname, sinfo in config["strains"].items():
    _STRAIN_BY_PREFIX[sinfo["prefix"]] = {
        "strain": sname,
        "accession": sinfo["accession"],
        "chrom_len": sinfo.get("chrom_len"),
    }

# Longest prefix first so "CJB111" wins over any shorter accidental match.
_PREFIXES_SORTED = sorted(_STRAIN_BY_PREFIX, key=len, reverse=True)

IP_TAG = config["naming"]["ip_tag"]
CTRL_TAG = config["naming"]["control_tag"]
BAM_SUF = config["naming"]["bam_suffix"]
R1_SUF = config["naming"]["r1_suffix"]
R2_SUF = config["naming"]["r2_suffix"]


def _strain_for_sample_token(token: str):
    """Given e.g. 'A9091' or 'CJB111C', return strain-info dict + replicate id."""
    for pref in _PREFIXES_SORTED:
        if token.startswith(pref):
            rep = token[len(pref):]              # "1","2","3" or "C"
            info = dict(_STRAIN_BY_PREFIX[pref])
            info["replicate"] = rep
            return info
    raise ValueError(
        f"Sample token '{token}' matches no strain prefix in config['strains']. "
        f"Known prefixes: {_PREFIXES_SORTED}"
    )


def _parse_bam_name(relpath: str):
    """Parse one BAM archive path into a sample metadata record."""
    base = PurePosixPath(relpath).name
    prefix = base[: -len(BAM_SUF)]               # strip _FilteredBAM.bam
    # {coreID}_{tag}_{antibody}_{sample}_-_-_-_{assay}
    fields = prefix.split("_")
    core_id, tag, antibody, sample_tok = fields[0], fields[1], fields[2], fields[3]
    sinfo = _strain_for_sample_token(sample_tok)
    is_ip = (tag == IP_TAG)
    return {
        "id": prefix,                            # unique wildcard token
        "core_id": core_id,
        "tag": tag,
        "antibody": antibody,
        "sample_token": sample_tok,
        "strain": sinfo["strain"],
        "accession": sinfo["accession"],
        "chrom_len": sinfo["chrom_len"],
        "replicate": sinfo["replicate"],
        "kind": "IP" if is_ip else "control",
        # archive-relative paths (same stem for bam / fastq)
        "arc_bam": relpath,
        "arc_r1": relpath[: -len(BAM_SUF)] + R1_SUF,
        "arc_r2": relpath[: -len(BAM_SUF)] + R2_SUF,
    }


def _build_sample_table(zip_path: str):
    if not os.path.exists(zip_path):
        raise FileNotFoundError(
            f"input_zip not found: {zip_path}\n"
            f"Set config['input_zip'] to the sequencing-core .zip."
        )
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
    bams = [n for n in names if n.endswith(BAM_SUF)]
    if not bams:
        raise ValueError(f"No '*{BAM_SUF}' entries found in {zip_path}.")
    table = {}
    for b in sorted(bams):
        rec = _parse_bam_name(b)
        table[rec["id"]] = rec
    return table


# ---------------------------------------------------------------------------
# Entry point: archive or FASTQ
#
# Discovery must produce IDENTICAL records either way, because every rule
# downstream addresses samples by the CEGR basename and reads strain, kind and
# replicate out of these fields. The only differences are where the BAM comes
# from (extracted vs aligned) and where the reads live (inside the archive vs
# supplied by the user).
# ---------------------------------------------------------------------------
def _build_sample_table_from_sheet(sheet_path: str, fastq_dir: str):
    """Discover samples from config/samples_fastq.tsv.

    sample_id must reproduce the CEGR basename exactly, because it is parsed
    with the same splitter used for archive entries. That keeps one convention
    rather than two, and means a typo in the sheet fails loudly here instead of
    producing a silently mislabelled sample.
    """
    if not os.path.exists(sheet_path):
        raise FileNotFoundError(
            f"input_mode is 'fastq' but the sample sheet was not found: "
            f"{sheet_path}\nCopy config/samples_fastq.tsv.example and fill it in."
        )
    table = {}
    with open(sheet_path) as fh:
        # Skip leading comments and blank lines: the shipped example carries an
        # explanatory header block, and a user editing it will keep those lines.
        hdr = None
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            hdr = ln.rstrip("\n").split("\t")
            break
        if hdr is None:
            raise ValueError(f"{sheet_path} contains no header row.")
        need = {"sample_id", "strain", "kind", "fastq_r1", "fastq_r2"}
        missing = need - set(hdr)
        if missing:
            raise ValueError(f"{sheet_path} lacks column(s): {sorted(missing)}")
        idx = {c: n for n, c in enumerate(hdr)}
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            sid = f[idx["sample_id"]].strip()

            # Same parse as _parse_bam_name, minus the suffix stripping.
            fields = sid.split("_")
            if len(fields) < 4:
                raise ValueError(
                    f"sample_id '{sid}' does not follow the CEGR convention "
                    f"{{coreID}}_{{tag}}_{{antibody}}_{{sample}}_-_-_-_{{assay}}"
                )
            core_id, tag, antibody, sample_tok = fields[0], fields[1], fields[2], fields[3]
            sinfo = _strain_for_sample_token(sample_tok)

            # Cross-check the sheet against what the name implies. A mismatch
            # means the sheet and the naming disagree, which would mislabel a
            # sample everywhere downstream.
            declared = f[idx["strain"]].strip()
            if declared != sinfo["strain"]:
                raise ValueError(
                    f"{sid}: sheet says strain '{declared}' but the sample "
                    f"token '{sample_tok}' implies '{sinfo['strain']}'"
                )
            kind = f[idx["kind"]].strip()
            if kind not in ("IP", "control"):
                raise ValueError(f"{sid}: kind must be 'IP' or 'control', got '{kind}'")
            implied = "IP" if tag == IP_TAG else "control"
            if kind != implied:
                raise ValueError(
                    f"{sid}: sheet says kind '{kind}' but tag '{tag}' implies "
                    f"'{implied}'"
                )

            table[sid] = {
                "id": sid,
                "core_id": core_id,
                "tag": tag,
                "antibody": antibody,
                "sample_token": sample_tok,
                "strain": sinfo["strain"],
                "accession": sinfo["accession"],
                "chrom_len": sinfo["chrom_len"],
                "replicate": sinfo["replicate"],
                "kind": kind,
                # BAMs are produced by align.smk directly under RAW, so the
                # "archive-relative" path is just the basename.
                "arc_bam": sid + BAM_SUF,
                # Reads are user-supplied and live outside RAW, so these are
                # full paths and raw_r1/raw_r2 return them unmodified.
                "arc_r1": None,
                "arc_r2": None,
                "fq_r1": os.path.join(fastq_dir, f[idx["fastq_r1"]].strip()),
                "fq_r2": os.path.join(fastq_dir, f[idx["fastq_r2"]].strip()),
            }
    if not table:
        raise ValueError(f"No sample rows found in {sheet_path}.")
    return table


INPUT_MODE = config.get("input_mode", "zip")
if INPUT_MODE == "fastq":
    _al = config.get("alignment", {})
    SAMPLES = _build_sample_table_from_sheet(
        _al.get("sample_sheet", "config/samples_fastq.tsv"),
        _al.get("fastq_dir", "data/fastq"),
    )
elif INPUT_MODE == "zip":
    SAMPLES = _build_sample_table(INPUT_ZIP)
else:
    raise ValueError(f"config['input_mode'] must be 'zip' or 'fastq', "
                     f"got '{INPUT_MODE}'")

SAMPLE_IDS = sorted(SAMPLES)
IP_IDS = [s for s in SAMPLE_IDS if SAMPLES[s]["kind"] == "IP"]
CTRL_IDS = [s for s in SAMPLE_IDS if SAMPLES[s]["kind"] == "control"]
STRAINS = sorted({SAMPLES[s]["strain"] for s in SAMPLE_IDS})

# strain -> ordered sample ids (IPs then control) for grouped QC (fingerprint)
def strain_ids(strain, kind=None):
    ids = [s for s in SAMPLE_IDS if SAMPLES[s]["strain"] == strain]
    if kind:
        ids = [s for s in ids if SAMPLES[s]["kind"] == kind]
    return sorted(ids)


# --- path helpers (post-extraction locations) --------------------------------
def raw_bam(sample):
    return os.path.join(RAW, SAMPLES[sample]["arc_bam"])

def raw_r1(sample):
    # In fastq mode the reads are supplied by the user and are not extracted
    # into RAW, so the sheet path is returned as-is.
    if INPUT_MODE == "fastq":
        return SAMPLES[sample]["fq_r1"]
    return os.path.join(RAW, SAMPLES[sample]["arc_r1"])

def raw_r2(sample):
    if INPUT_MODE == "fastq":
        return SAMPLES[sample]["fq_r2"]
    return os.path.join(RAW, SAMPLES[sample]["arc_r2"])

def raw_bai(sample):
    # BAM index lives adjacent to the BAM (required by deepTools/ChExMix)
    return raw_bam(sample) + ".bai"


# Constrain wildcards to the exact set of discovered sample ids so the
# underscores/dashes in the names never cause ambiguous matches.
wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLE_IDS),
    strain="|".join(re.escape(s) for s in STRAINS),


# All raw files that stage_delivery must produce (extraction targets).
def all_raw_targets():
    """Files that must exist before QC can run.

    In zip mode all three are produced by extraction. In fastq mode the BAM is
    produced by align.smk and the reads are inputs the user has already
    downloaded, so only the BAM is a build target.
    """
    t = []
    for s in SAMPLE_IDS:
        t.append(raw_bam(s))
        if INPUT_MODE == "zip":
            t += [raw_r1(s), raw_r2(s)]
    return t


# =============================================================================
# Helpers for the PEAKS stage
# =============================================================================
REF = config["reference"]["outdir"]
PEAKS = os.path.join(RESULTS, "03_peaks")
# Peak calls are segregated by CALLER so that a second caller can be added
# without either overwriting the other. PEAKS remains the parent directory.
PEAKS_CHEXMIX = os.path.join(PEAKS, "chexmix")
PEAKS_GEM     = os.path.join(PEAKS, "gem")


def strain_accession(strain):
    """The single reference accession for a strain (asserts consistency)."""
    accs = {SAMPLES[s]["accession"] for s in strain_ids(strain)}
    assert len(accs) == 1, f"{strain} maps to multiple accessions: {accs}"
    return accs.pop()

def strain_assembly(strain):
    return config["strains"][strain]["assembly"]

def all_assemblies():
    return sorted({v["assembly"] for v in config["strains"].values()})

# reference files are stored per ASSEMBLY accession
def ref_fna(asm):     return os.path.join(REF, asm, f"{asm}.fna")
def ref_gff(asm):     return os.path.join(REF, asm, f"{asm}.gff")
def ref_geninfo(asm): return os.path.join(REF, asm, f"{asm}.geninfo")
def ref_seqdir(asm):  return os.path.join(REF, asm, "seq")