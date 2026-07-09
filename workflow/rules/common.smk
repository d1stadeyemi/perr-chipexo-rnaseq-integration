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
INPUT_ZIP = config["input_zip"]

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


SAMPLES = _build_sample_table(INPUT_ZIP)
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
    return os.path.join(RAW, SAMPLES[sample]["arc_r1"])

def raw_r2(sample):
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
    t = []
    for s in SAMPLE_IDS:
        t += [raw_bam(s), raw_r1(s), raw_r2(s)]
    return t


# =============================================================================
# Helpers for the PEAKS stage
# =============================================================================
REF = config["reference"]["outdir"]
PEAKS = os.path.join(RESULTS, "03_peaks")


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