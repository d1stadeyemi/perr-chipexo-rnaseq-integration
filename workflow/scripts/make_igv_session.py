"""make_igv_session.py

Generate IGV session files and a batch script.

WHY THIS IS WORTH SCRIPTING
  An IGV session is XML, so every display setting that affects what a reader
  concludes can be written by a script rather than clicked by a person:

    autoScale = false      With autoscale ON, a weak CJB111 signal and a strong
                           A909 signal render IDENTICALLY, because each is scaled
                           to its own maximum. This is the most common way a
                           browser figure misleads, and a reviewer who uses IGV
                           will spot it. Explicit, shared DataRange instead.

    windowFunction = max   IGV's default is Mean, which averages away
                           single-base spikes -- and single-base spikes ARE the
                           ChIP-exo signal. At the crosslink this is the
                           difference between seeing a footprint and not.

    fixed track order      Forward, reverse, control, peaks, motifs, genes.
                           The control is never omitted: a peak shown without a
                           flat control beneath it is not evidence.

  Opening one file therefore reproduces a fully configured browser, and the
  settings live in version control rather than in someone's memory.

TWO SESSIONS PER STRAIN x DEPTH
  <strain>_<depth>.xml           forward and reverse stacked, both positive
  <strain>_<depth>_mirrored.xml  reverse negated and drawn below the axis, which
                                 is the canonical ChIP-exo presentation and makes
                                 the strand offset read at a glance

A NOTE ON THE BATCH SCRIPT
  IGV batch scripts cannot set data ranges. The batch file therefore LOADS THE
  SESSION first, so the ranges come from the XML, and then only navigates and
  snapshots. Loading tracks individually in a batch script would silently lose
  every display setting.
"""
import os
from xml.sax.saxutils import escape

I = snakemake.input                          # noqa: F821
O = snakemake.output                         # noqa: F821
strain = snakemake.params.strain             # noqa: F821
depth = snakemake.params.depth               # noqa: F821
cfg = dict(snakemake.params.cfg)             # noqa: F821
logp = snakemake.log[0]                      # noqa: F821

os.makedirs(os.path.dirname(O.session), exist_ok=True)
os.makedirs(os.path.dirname(O.batch), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

A = os.path.abspath
COL = cfg.get("colors", {})
R_IP = float(cfg.get("range_ip", 150))
R_CT = float(cfg.get("range_ctrl", 150))
DEPTH_LABEL = "full depth" if depth == "full_depth" else "DEPTH-MATCHED"


def data_track(path, name, colour, dmin, dmax, height=70):
    """A bigWig track with an explicit range and autoscale disabled."""
    return f"""        <Track attributeKey="{escape(name)}" autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" color="{colour}" fontSize="10" height="{height}" id="{escape(A(path))}" name="{escape(name)}" renderer="BAR_CHART" visible="true" windowFunction="max">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="{dmax}" minimum="{dmin}" type="LINEAR"/>
        </Track>
"""


def feature_track(path, name, colour, height=28, display="COLLAPSED"):
    return f"""        <Track attributeKey="{escape(name)}" clazz="org.broad.igv.track.FeatureTrack" color="{colour}" displayMode="{display}" fontSize="10" height="{height}" id="{escape(A(path))}" name="{escape(name)}" visible="true"/>
"""


def build(mirrored, locus):
    """Assemble one session XML."""
    resources, data, feats = [], [], []

    def res(p):
        resources.append(f'        <Resource path="{escape(A(p))}"/>\n')

    # --- coverage -----------------------------------------------------------
    res(I.fwd)
    data.append(data_track(I.fwd, f"{strain} IP forward 5\u2032 ({DEPTH_LABEL})",
                           COL.get("forward", "0,90,160"), 0.0, R_IP))
    if mirrored:
        res(I.neg)
        data.append(data_track(I.neg, f"{strain} IP reverse 5\u2032 (mirrored)",
                               COL.get("reverse", "168,71,28"), -R_IP, 0.0))
    else:
        res(I.rev)
        data.append(data_track(I.rev, f"{strain} IP reverse 5\u2032 ({DEPTH_LABEL})",
                               COL.get("reverse", "168,71,28"), 0.0, R_IP))
    res(I.ctrl)
    data.append(data_track(I.ctrl, f"{strain} no-tag control ({DEPTH_LABEL})",
                           COL.get("control", "150,158,165"), 0.0, R_CT, height=50))

    # --- features -----------------------------------------------------------
    for path, name, key in (
        (I.peaks, f"{strain} ChExMix peaks (full depth)", "peaks"),
        (I.peaks_ds, f"{strain} peaks surviving depth match", "peaks"),
        (I.motif_gbs, f"{strain} GBS PerR motif (FIMO p<1e-4)", "motif_gbs"),
        (I.motif_gas, f"{strain} GAS Per box (reference only)", "motif_gas"),
        (I.regulon, f"{strain} direct regulon", "peaks"),
    ):
        res(path)
        feats.append(feature_track(path, name, COL.get(key, "31,41,51")))

    res(I.gff)
    feats.append(feature_track(I.gff, f"{strain} genes (RefSeq)", "60,70,80", height=40))

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Auto-generated by workflow/rules/igv.smk. Do not edit by hand: it is
     regenerated on every `snakemake igv` run. Tweak in IGV, then File > Save
     Session to a NEW filename if you want to keep changes.

     depth  = {depth}   ({'browse / single-strain figures' if depth == 'full_depth' else 'CROSS-STRAIN claims only'})
     ranges = IP 0-{R_IP:g}, control 0-{R_CT:g}, autoscale OFF, windowFunction = max
-->
<Session genome="{escape(A(I.fna))}" hasGeneTrack="false" hasSequenceTrack="true" locus="{escape(locus)}" version="8">
    <Resources>
{''.join(resources)}    </Resources>
    <Panel height="520" name="DataPanel" width="1200">
{''.join(data)}    </Panel>
    <Panel height="260" name="FeaturePanel" width="1200">
{''.join(feats)}    </Panel>
    <PanelLayout dividerFractions="0.66"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
"""


# ---------------------------------------------------------------------------
# opening locus: the first flagship panel for this strain, else the whole chrom
# ---------------------------------------------------------------------------
def first_locus():
    try:
        with open(I.flagship) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) > idx["igv_context"] and f[idx["strain"]] == strain \
                        and f[idx["igv_context"]]:
                    return f[idx["igv_context"]]
    except Exception:  # noqa: BLE001
        pass
    return "All"


LOCUS = first_locus()

with open(O.session, "w") as fh:
    fh.write(build(False, LOCUS))
with open(O.session_mirror, "w") as fh:
    fh.write(build(True, LOCUS))


# ---------------------------------------------------------------------------
# batch script: load the session, then navigate and snapshot
# ---------------------------------------------------------------------------
panels = []
try:
    with open(I.flagship) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr) or f[idx["strain"]] != strain:
                continue
            label = f[idx["panel"]].replace(" ", "_").replace("/", "-")
            label = "".join(ch for ch in label if ch.isalnum() or ch in "._-")
            if f[idx["igv_context"]]:
                panels.append((f[idx["igv_context"]], f"{f[idx['order']]}_{label}_context"))
            if f[idx["igv_zoom"]]:
                panels.append((f[idx["igv_zoom"]], f"{f[idx['order']]}_{label}_zoom"))
except Exception:  # noqa: BLE001
    pass

snapdir = A(os.path.join(os.path.dirname(os.path.dirname(O.session)), "panels", f"{strain}_{depth}"))
lines = [
    "# Auto-generated IGV batch script.",
    "# Run:  igv.sh --batch <this file>        (headless needs xvfb-run)",
    "#",
    "# The SESSION is loaded first so that data ranges, colours and the window",
    "# function come from the XML -- batch scripts cannot set those, and loading",
    "# tracks individually here would silently lose every display setting.",
    "new",
    f"load {A(O.session)}",
    "maxPanelHeight 900",
    f"snapshotDirectory {snapdir}",
]
for locus, name in panels:
    lines += [f"goto {locus}", f"snapshot {strain}_{depth}_{name}.png"]
lines += ["exit", ""]

os.makedirs(snapdir, exist_ok=True)
with open(O.batch, "w") as fh:
    fh.write("\n".join(lines))

with open(logp, "w") as lg:
    lg.write(f"strain={strain} depth={depth}\n")
    lg.write(f"opening locus: {LOCUS}\n")
    lg.write(f"sessions: {O.session}\n          {O.session_mirror}\n")
    lg.write(f"batch: {O.batch}  ({len(panels)} snapshots)\n")
    lg.write(f"IP range 0-{R_IP:g}, control 0-{R_CT:g}, autoscale OFF, windowFunction=max\n")
    if not panels:
        lg.write("WARNING: no flagship panels matched this strain; the batch "
                 "script will only open the session.\n")
