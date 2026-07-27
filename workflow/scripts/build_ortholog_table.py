"""Reciprocal best hits (RBH) between two proteomes -> locus_tag ortholog pairs.

Why RBH: with two closely related strains of the same species, RBH is the
standard, conservative definition of one-to-one orthology. It is symmetric by
construction, so the A909->CJB111 and CJB111->A909 regulon comparisons agree.
RBH is an OPERATIONAL proxy for orthology, not a phylogenetic one; its known
failure mode is recent within-genome duplication, where two copies compete for
the reciprocal slot and one true ortholog can be dropped. At GBS intra-species
distances that is rare.

Protein accessions are mapped back to locus tags using the CDS features of each
GFF (which carry both protein_id= and locus_tag=). NOTE: identical proteins in
a genome share one WP_ accession, so a protein may map to several locus tags;
such ambiguous pairs are flagged rather than silently dropped. This is why the
number of OUTPUT ROWS exceeds the number of RBH pairs -- one protein-level hit
can expand into several locus-tag rows (typically IS elements / transposases).

QCOVERAGE > 100 IS NOT A BUG. Coverage is alignment_length / query_length, and
BLAST's `length` counts gap columns, so a gapped alignment can span more
columns than the query has residues. Values slightly above 100 are expected.

LOW-IDENTITY FLAG (v2): between two strains of one species, genuine orthologs
sit above ~95% identity. A reciprocal best hit far below that is a candidate
PARALOG CAPTURE or a mosaic/recombinant locus (e.g. variable surface antigens
such as the C5a peptidase scpA/B), and must be reviewed individually before it
enters a figure. Such rows are flagged LOW_IDENTITY rather than dropped, so the
judgement stays with the analyst and is visible in the output.
"""
import os
import re

fwd = snakemake.input.fwd                # noqa: F821  A -> B blastp (outfmt 6)
rev = snakemake.input.rev                # noqa: F821  B -> A blastp
gff_a = snakemake.input.gff_a            # noqa: F821
gff_b = snakemake.input.gff_b            # noqa: F821
out_tsv = snakemake.output.tsv           # noqa: F821
sa = snakemake.params.strain_a           # noqa: F821
sb = snakemake.params.strain_b           # noqa: F821
min_id = float(snakemake.params.min_identity)     # noqa: F821
min_cov = float(snakemake.params.min_coverage)    # noqa: F821
flag_id = float(snakemake.params.flag_identity)   # noqa: F821
logp = snakemake.log[0]                  # noqa: F821

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
os.makedirs(os.path.dirname(logp), exist_ok=True)


def prot2locus(gff):
    """protein_id -> [locus_tag, ...] from CDS features."""
    m = {}
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            pid = re.search(r"protein_id=([^;]+)", f[8])
            lt = re.search(r"locus_tag=([^;]+)", f[8])
            if pid and lt:
                m.setdefault(pid.group(1), [])
                if lt.group(1) not in m[pid.group(1)]:
                    m[pid.group(1)].append(lt.group(1))
    return m


def best_hits(path):
    """query -> (subject, pident, qcov, bitscore), top-scoring line per query.
    blastp outfmt: qseqid sseqid pident length qlen slen evalue bitscore
    Bitscore (not raw score or e-value) is the ranking key: it is comparable
    across searches against databases of different size."""
    best = {}
    with open(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            pid, alen, qlen = float(f[2]), int(f[3]), int(f[4])
            bits = float(f[7])
            cov = 100.0 * alen / qlen if qlen else 0.0
            if pid < min_id or cov < min_cov:
                continue
            if q not in best or bits > best[q][3]:
                best[q] = (s, pid, cov, bits)
    return best


pa, pb = prot2locus(gff_a), prot2locus(gff_b)
fb, rb = best_hits(fwd), best_hits(rev)

rows, n_rbh, n_ambig, n_lowid = [], 0, 0, 0
for q, (s, pid, cov, _) in fb.items():
    if s in rb and rb[s][0] == q:                 # reciprocal
        n_rbh += 1
        la, lb = pa.get(q, []), pb.get(s, [])
        if not la or not lb:
            continue
        ambiguous = (len(la) > 1 or len(lb) > 1)
        if ambiguous:
            n_ambig += 1
        flag = "LOW_IDENTITY" if pid < flag_id else ""
        if flag:
            n_lowid += 1
        for x in la:
            for y in lb:
                rows.append([x, y, q, s, f"{pid:.1f}", f"{cov:.1f}",
                             "ambiguous" if ambiguous else "1:1", flag])

with open(out_tsv, "w") as ob:
    ob.write(f"{sa}_locus_tag\t{sb}_locus_tag\t{sa}_protein\t{sb}_protein\t"
             f"pident\tqcoverage\tmapping\tflag\n")
    for r in sorted(rows):
        ob.write("\t".join(r) + "\n")

with open(logp, "w") as lg:
    lg.write(f"{sa} proteins with locus tags: {len(pa)}\n")
    lg.write(f"{sb} proteins with locus tags: {len(pb)}\n")
    lg.write(f"reciprocal best hits: {n_rbh}  (ambiguous protein->locus: {n_ambig})\n")
    lg.write(f"ortholog rows written: {len(rows)}\n")
    lg.write(f"RBH pairs flagged LOW_IDENTITY (<{flag_id}%): {n_lowid}"
             f"  -- review individually before use\n")
    lg.write(f"filters: pident>={min_id} qcoverage>={min_cov}\n")