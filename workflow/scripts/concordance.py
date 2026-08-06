"""concordance.py

Cross-caller agreement between ChExMix and GEM.

WHY SEVERAL TOLERANCES RATHER THAN ONE
  The two callers segment dense regions differently. GEM resolved the A909 perR
  locus into ~5 components across 450 bp where ChExMix reported 2, so at a 25 bp
  tolerance three GEM events look "unique" when they are in fact inside a region
  ChExMix also called. Reporting a single number would therefore either
  understate agreement (strict) or overstate it (loose). Both the event-level
  and the region-level view are emitted so the reader can see the difference.

REGION-LEVEL MATCHING
  A GEM event is counted as region-matched if it falls inside the span of a
  ChExMix peak CLUSTER -- consecutive ChExMix peaks separated by less than
  `cluster_gap`. This is the number that answers "did the two callers find the
  same loci", as distinct from "did they place events identically".

GREEDY ONE-TO-ONE PAIRING
  Event-level matching pairs each GEM event with its nearest unused ChExMix
  peak. Without the one-to-one constraint, five GEM components inside one
  ChExMix peak would all count as matches and inflate the figure.
"""
import os

chexmix_path = snakemake.input.chexmix          # noqa: F821
gem_path = snakemake.input.gem                  # noqa: F821
out = snakemake.output.tsv                      # noqa: F821
detail = snakemake.output.detail                # noqa: F821
strain = snakemake.params.strain                # noqa: F821
TOLS = [int(t) for t in snakemake.params.tolerances]   # noqa: F821
FLAG = dict(snakemake.params.flagship or {})    # noqa: F821
logp = snakemake.log[0]                         # noqa: F821

os.makedirs(os.path.dirname(out), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)

CLUSTER_GAP = 500


def read_events(path):
    """-> [(pos, score)] from either a ChExMix .events or a GEM GPS/GEM_events
    file. Both put 'chrom:pos' in column 1; the score column differs, so it is
    read defensively and only used for reporting."""
    rows = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(("#", "Position")) or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            try:
                pos = int(f[0].split(":")[1])
            except (IndexError, ValueError):
                continue
            sc = ""
            for i in (3, 5):
                if len(f) > i:
                    try:
                        sc = float(f[i].strip())
                        break
                    except ValueError:
                        pass
            rows.append((pos, sc))
    return sorted(rows)


CX = read_events(chexmix_path)
GM = read_events(gem_path)
cx_pos = [p for p, _ in CX]


def clusters(positions, gap):
    """Merge positions closer than `gap` into spans."""
    if not positions:
        return []
    out, s, e = [], positions[0], positions[0]
    for p in positions[1:]:
        if p - e <= gap:
            e = p
        else:
            out.append((s, e))
            s = e = p
    out.append((s, e))
    return out


CX_CL = clusters(cx_pos, CLUSTER_GAP)


def pair(a, b, tol):
    """Greedy one-to-one nearest pairing."""
    used, pairs = set(), []
    for pa, _ in a:
        best = None
        for i, (pb, _) in enumerate(b):
            if i in used or abs(pa - pb) > tol:
                continue
            if best is None or abs(pa - pb) < abs(pa - b[best][0]):
                best = i
        if best is not None:
            used.add(best)
            pairs.append((pa, b[best][0], abs(pa - b[best][0])))
    return pairs


rows = []
for tol in TOLS:
    g2c = pair(GM, CX, tol)
    c2g = pair(CX, GM, tol)
    rows.append(dict(
        tolerance_bp=tol,
        gem_events=len(GM),
        chexmix_events=len(CX),
        gem_matched=len(g2c),
        gem_matched_pct=round(100.0 * len(g2c) / len(GM), 1) if GM else "",
        chexmix_matched=len(c2g),
        chexmix_matched_pct=round(100.0 * len(c2g) / len(CX), 1) if CX else "",
        median_offset_bp=(sorted(d for _, _, d in g2c)[len(g2c) // 2] if g2c else ""),
    ))

# region-level: is each GEM event inside a ChExMix peak cluster?
in_region = sum(1 for p, _ in GM if any(s - 100 <= p <= e + 100 for s, e in CX_CL))

with open(out, "w") as ob:
    cols = list(rows[0].keys())
    ob.write("strain\t" + "\t".join(cols) + "\n")
    for r in rows:
        ob.write(strain + "\t" + "\t".join(str(r[c]) for c in cols) + "\n")
    ob.write(f"{strain}\tREGION_LEVEL\t{len(GM)}\t{len(CX)}\t{in_region}\t"
             f"{round(100.0*in_region/len(GM),1) if GM else ''}\t\t\t\n")

with open(detail, "w") as db:
    db.write("gem_pos\tnearest_chexmix\toffset_bp\tin_chexmix_cluster\n")
    for p, _ in GM:
        near = min(cx_pos, key=lambda q: abs(p - q)) if cx_pos else ""
        off = abs(p - near) if near != "" else ""
        inc = any(s - 100 <= p <= e + 100 for s, e in CX_CL)
        db.write(f"{p}\t{near}\t{off}\t{inc}\n")

with open(logp, "w") as lg:
    lg.write(f"strain={strain}\n")
    lg.write(f"ChExMix events: {len(CX)}   GEM events: {len(GM)}\n")
    lg.write(f"ChExMix peak clusters (gap<{CLUSTER_GAP}bp): {len(CX_CL)}\n")
    lg.write(f"GEM events inside a ChExMix cluster: {in_region}/{len(GM)}\n\n")
    for r in rows:
        lg.write(f"  tol {r['tolerance_bp']:4d} bp: GEM {r['gem_matched']}/{len(GM)} "
                 f"({r['gem_matched_pct']}%)  median offset {r['median_offset_bp']} bp\n")
    if FLAG:
        lg.write("\nflagship loci:\n")
        for name, pos in FLAG.items():
            hit = min((abs(int(pos) - p), p) for p, _ in GM) if GM else (None, None)
            lg.write(f"  {name:22s} {pos}: "
                     + (f"GEM {hit[1]} ({hit[0]} bp away)" if hit[0] is not None
                        and hit[0] <= 200 else "NOT recovered by GEM") + "\n")