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
                           single-base spikes, and single-base spikes ARE the
                           ChIP-exo signal. At the crosslink this is the
                           difference between seeing a footprint and not.

    fixed track order      Forward, reverse, control, peaks, motifs, genes.
                           The control is never omitted: a peak shown without a
                           flat control beneath it is not evidence.

-----------------------------------------------------------------------------
v2  PATH PORTABILITY
-----------------------------------------------------------------------------
  v1 wrote ABSOLUTE paths (os.path.abspath). Those resolve only on the machine
  that generated them, so copying results/10_igv/ to a laptop produced a session
  in which every track loaded SILENTLY EMPTY. IGV does not warn about a missing
  resource; it simply shows nothing.

  IGV resolves RELATIVE resource paths against the SESSION FILE's own directory,
  so paths are now written relative to that directory and the tree can be moved
  anywhere. The reference FASTA/FAI/GFF are staged into results/10_igv/genome/
  by rule igv_stage_genome, so that

      tar czf igv_bundle.tar.gz results/10_igv

  is a complete bundle with nothing outside it to chase.

  Set igv.path_style: "absolute" in config to restore the old behaviour.

-----------------------------------------------------------------------------
v3  XML COMMENT SAFETY  (see xml_comment below)
-----------------------------------------------------------------------------
  v2 emitted a header comment containing "--", which XML 1.0 section 2.5
  forbids inside comment content. Every session failed to load with

      org.xml.sax.SAXParseException; The string "--" is not permitted
      within comments.

  The header is assembled by interpolating variables (depth, path style, ranges,
  notes), so hand-fixing the offending literals would not prevent recurrence:
  any future value containing a double hyphen would break all sessions again,
  with the same opaque parser error. The text is therefore passed through
  xml_comment(), which guarantees a valid comment regardless of input.

A NOTE ON THE BATCH SCRIPT
  IGV batch scripts cannot set data ranges. The batch file therefore LOADS THE
  SESSION first, so ranges come from the XML, and then only navigates and
  snapshots. Loading tracks individually would silently lose every setting.
  Batch paths resolve against the WORKING DIRECTORY, not the batch file, so the
  header states which directory to run it from.
"""
import os
import re
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

COL = cfg.get("colors", {})
R_IP = float(cfg.get("range_ip", 150))
R_CT = float(cfg.get("range_ctrl", 150))
STYLE = str(cfg.get("path_style", "relative")).lower()
DEPTH_LABEL = "full depth" if depth == "full_depth" else "DEPTH-MATCHED"

# Both sessions (plain and mirrored) live in the same directory, so one anchor
# serves both.
SESSION_DIR = os.path.dirname(os.path.abspath(O.session))
REPO_ROOT = os.getcwd()


def xml_comment(text):
    """Return text that is guaranteed safe inside an XML comment.

    XML 1.0 section 2.5 imposes two constraints that are easy to violate by
    accident when a comment is built from interpolated values:

      1. the two-character sequence "--" may not appear in comment content
      2. a comment may not end on "-" (because "--->" would be ambiguous)

    Runs of hyphens are collapsed to an en dash rather than deleted, so the text
    still reads as intended: "RELATIVE -- this tree" becomes
    "RELATIVE \u2013 this tree" instead of losing the separator entirely. The en
    dash is ordinary UTF-8 and the document declares that encoding.

    This is deliberately defensive. Sanitising here means no configuration value
    or generated path can ever produce an unparseable session, which is a better
    guarantee than remembering not to type a double hyphen."""
    out = re.sub(r"-{2,}", "\u2013", text)
    return out.rstrip().rstrip("-")


def P(path):
    """Path as written into the session XML.

    relative -> resolved by IGV against the session file's directory, so the
                tree can be copied anywhere
    absolute -> resolves only on this machine
    """
    ap = os.path.abspath(path)
    if STYLE == "absolute":
        return ap
    # posix separators: IGV accepts them on every platform
    return os.path.relpath(ap, SESSION_DIR).replace(os.sep, "/")


def B(path):
    """Path as written into the BATCH script. Batch paths resolve against the
    working directory, so these are relative to the repository root."""
    ap = os.path.abspath(path)
    if STYLE == "absolute":
        return ap
    return os.path.relpath(ap, REPO_ROOT).replace(os.sep, "/")


def data_track(path, name, colour, dmin, dmax, height=70):
    """A bigWig track with an explicit range and autoscale disabled."""
    return (
        '        <Track attributeKey="%s" autoScale="false"'
        ' clazz="org.broad.igv.track.DataSourceTrack" color="%s" fontSize="10"'
        ' height="%d" id="%s" name="%s" renderer="BAR_CHART" visible="true"'
        ' windowFunction="max">\n'
        '            <DataRange baseline="0.0" drawBaseline="true"'
        ' flipAxis="false" maximum="%s" minimum="%s" type="LINEAR"/>\n'
        '        </Track>\n'
    ) % (escape(name), colour, height, escape(P(path)), escape(name), dmax, dmin)


def feature_track(path, name, colour, height=28, display="COLLAPSED"):
    return (
        '        <Track attributeKey="%s"'
        ' clazz="org.broad.igv.track.FeatureTrack" color="%s" displayMode="%s"'
        ' fontSize="10" height="%d" id="%s" name="%s" visible="true"/>\n'
    ) % (escape(name), colour, display, height, escape(P(path)), escape(name))


def build(mirrored, locus):
    """Assemble one session XML."""
    resources, data, feats = [], [], []

    def res(p):
        resources.append('        <Resource path="%s"/>\n' % escape(P(p)))

    # --- coverage -----------------------------------------------------------
    res(I.fwd)
    data.append(data_track(I.fwd, "%s IP forward 5\u2032 (%s)" % (strain, DEPTH_LABEL),
                           COL.get("forward", "0,90,160"), 0.0, R_IP))
    if mirrored:
        res(I.neg)
        data.append(data_track(I.neg, "%s IP reverse 5\u2032 (mirrored)" % strain,
                               COL.get("reverse", "168,71,28"), -R_IP, 0.0))
    else:
        res(I.rev)
        data.append(data_track(I.rev, "%s IP reverse 5\u2032 (%s)" % (strain, DEPTH_LABEL),
                               COL.get("reverse", "168,71,28"), 0.0, R_IP))
    res(I.ctrl)
    data.append(data_track(I.ctrl, "%s no-tag control (%s)" % (strain, DEPTH_LABEL),
                           COL.get("control", "150,158,165"), 0.0, R_CT, height=50))

    # --- features -----------------------------------------------------------
    for path, name, key in (
        (I.peaks, "%s ChExMix peaks (full depth)" % strain, "peaks"),
        (I.peaks_ds, "%s peaks surviving depth match" % strain, "peaks"),
        (I.motif_gbs, "%s GBS PerR motif (FIMO p<1e-4)" % strain, "motif_gbs"),
        (I.motif_gas, "%s GAS Per box (reference only)" % strain, "motif_gas"),
        (I.regulon, "%s direct regulon" % strain, "peaks"),
    ):
        res(path)
        feats.append(feature_track(path, name, COL.get(key, "31,41,51")))

    res(I.gff)
    feats.append(feature_track(I.gff, "%s genes (RefSeq)" % strain, "60,70,80", height=40))

    # --- header comment, sanitised ------------------------------------------
    style_note = ("relative to this file's directory, so this tree can be copied anywhere"
                  if STYLE != "absolute" else
                  "absolute, so it resolves only on the machine that generated it")
    depth_note = ("browse / single-strain figures" if depth == "full_depth"
                  else "CROSS-STRAIN claims only")

    header = xml_comment(
        "Auto-generated by workflow/rules/igv.smk. Do not edit by hand: it is\n"
        "     regenerated on every `snakemake igv` run. Tweak in IGV, then\n"
        "     File > Save Session to a NEW filename if you want to keep changes.\n"
        "\n"
        "     depth  = %s   (%s)\n"
        "     ranges = IP 0 to %g, control 0 to %g, autoscale OFF, windowFunction = max\n"
        "     paths  = %s  (%s)\n"
        "\n"
        "     To move this to another machine:\n"
        "         tar czf igv_bundle.tar.gz results/10_igv\n"
        "     extract anywhere, then File > Open Session. Nothing outside the\n"
        "     bundle is referenced; the reference genome is staged into\n"
        "     results/10_igv/genome/."
        % (depth, depth_note, R_IP, R_CT, STYLE, style_note)
    )

    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        '<!-- %s -->\n'
        '<Session genome="%s" hasGeneTrack="false" hasSequenceTrack="true"'
        ' locus="%s" version="8">\n'
        '    <Resources>\n'
        '%s    </Resources>\n'
        '    <Panel height="520" name="DataPanel" width="1200">\n'
        '%s    </Panel>\n'
        '    <Panel height="260" name="FeaturePanel" width="1200">\n'
        '%s    </Panel>\n'
        '    <PanelLayout dividerFractions="0.66"/>\n'
        '    <HiddenAttributes>\n'
        '        <Attribute name="DATA FILE"/>\n'
        '        <Attribute name="DATA TYPE"/>\n'
        '        <Attribute name="NAME"/>\n'
        '    </HiddenAttributes>\n'
        '</Session>\n'
    ) % (header, escape(P(I.fna)), escape(locus),
         "".join(resources), "".join(data), "".join(feats))


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
# self-check: parse what we just wrote, so a malformed session can never leave
# this rule silently. IGV's error message is opaque; this one names the file.
# ---------------------------------------------------------------------------
import xml.dom.minidom  # noqa: E402

for _f in (O.session, O.session_mirror):
    try:
        xml.dom.minidom.parse(_f)
    except Exception as _e:  # noqa: BLE001
        raise SystemExit(
            "FATAL: generated session is not valid XML: %s\n  %s\n"
            "This should be impossible; xml_comment() is meant to prevent it. "
            "Please report the offending text." % (_f, _e))


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
                panels.append((f[idx["igv_context"]],
                               "%s_%s_context" % (f[idx["order"]], label)))
            if f[idx["igv_zoom"]]:
                panels.append((f[idx["igv_zoom"]],
                               "%s_%s_zoom" % (f[idx["order"]], label)))
except Exception:  # noqa: BLE001
    pass

snapdir_abs = os.path.abspath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(O.session))),
    "panels", "%s_%s" % (strain, depth)))
os.makedirs(snapdir_abs, exist_ok=True)

lines = [
    "# Auto-generated IGV batch script.",
    "#",
    "# PATH STYLE: %s" % STYLE,
    "# RUN FROM THE REPOSITORY ROOT:",
    "#   xvfb-run -a igv.sh --batch %s" % B(O.batch),
    "#",
    "# The SESSION is loaded first so that data ranges, colours and the window",
    "# function come from the XML. Batch scripts cannot set those, and loading",
    "# tracks individually here would silently lose every display setting.",
    "new",
    "load %s" % B(O.session),
    "maxPanelHeight 900",
    "snapshotDirectory %s" % B(snapdir_abs),
]
for locus, name in panels:
    lines += ["goto %s" % locus,
              "snapshot %s_%s_%s.png" % (strain, depth, name)]
lines += ["exit", ""]

with open(O.batch, "w") as fh:
    fh.write("\n".join(lines))

with open(logp, "w") as lg:
    lg.write("strain=%s depth=%s\n" % (strain, depth))
    lg.write("path_style=%s\n" % STYLE)
    lg.write("relative anchor (session dir): %s\n" % SESSION_DIR)
    lg.write("opening locus: %s\n" % LOCUS)
    lg.write("XML self-check: both sessions parsed OK\n")
    lg.write("--- paths as written into the XML ---\n")
    lg.write("  genome : %s\n" % P(I.fna))
    lg.write("  bigwig : %s\n" % P(I.fwd))
    lg.write("  peaks  : %s\n" % P(I.peaks))
    lg.write("sessions: %s\n          %s\n" % (O.session, O.session_mirror))
    lg.write("batch: %s  (%d snapshots)\n" % (O.batch, len(panels)))
    lg.write("IP range 0-%g, control 0-%g, autoscale OFF, windowFunction=max\n"
             % (R_IP, R_CT))
    if not panels:
        lg.write("WARNING: no flagship panels matched this strain; the batch "
                 "script will only open the session.\n")