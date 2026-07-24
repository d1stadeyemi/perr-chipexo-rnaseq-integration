"""Reciprocal best hits (RBH) between two proteomes -> locus_tag ortholog pairs.

Why RBH: with two closely related strains of the same species, RBH is the
standard, conservative definition of one-to-one orthology. It is symmetric by
construction, so the A909->CJB111 and CJB111->A909 regulon comparisons agree.

Protein accessions are mapped back to locus tags using the CDS features of each
GFF (which carry both protein_id= and locus_tag=). NOTE: identical proteins in
a genome share one WP_ accession, so a protein may map to several locus tags;
such ambiguous pairs are flagged rather than silently dropped.
"""
import os
import re
import csv

fwd     = snakemake.input.fwd            # noqa: F821  A -> B blastp (outfmt 6)
rev     = snakemake.input.rev            # noqa: F821  B -> A blastp
gff_a   = snakemake.input.gff_a          # noqa: F821
gff_b   = snakemake.input.gff_b          # noqa: F821
out_tsv = snakemake.output.tsv           # noqa: F821
sa      = snakemake.params.strain_a      # noqa: F821
sb      = snakemake.params.strain_b      # noqa: F821
min_id  = float(snakemake.params.min_identity)    # noqa: F821
min_cov = float(snakemake.params.min_coverage)    # noqa: F821
logp    = snakemake.log[0]               # noqa: F821

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
            lt  = re.search(r"locus_tag=([^;]+)", f[8])
            if pid and lt:
                m.setdefault(pid.group(1), [])
                if lt.group(1) not in m[pid.group(1)]:
                    m[pid.group(1)].append(lt.group(1))
    return m

def best_hits(path):
    """query -> (subject, pident, qcov) keeping the top-scoring line per query.
    blastp outfmt: qseqid sseqid pident length qlen slen evalue bitscore"""
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

rows, n_rbh, n_ambig = [], 0, 0
for q, (s, pid, cov, _) in fb.items():
    if s in rb and rb[s][0] == q:                 # reciprocal
        n_rbh += 1
        la, lb = pa.get(q, []), pb.get(s, [])
        if not la or not lb:
            continue
        if len(la) > 1 or len(lb) > 1:
            n_ambig += 1
        for x in la:
            for y in lb:
                rows.append([x, y, q, s, f"{pid:.1f}", f"{cov:.1f}",
                             "ambiguous" if (len(la) > 1 or len(lb) > 1) else "1:1"])

with open(out_tsv, "w") as ob:
    ob.write(f"{sa}_locus_tag\t{sb}_locus_tag\t{sa}_protein\t{sb}_protein\t"
             f"pident\tqcoverage\tmapping\n")
    for r in sorted(rows):
        ob.write("\t".join(r) + "\n")

with open(logp, "w") as lg:
    lg.write(f"{sa} proteins with locus tags: {len(pa)}\n")
    lg.write(f"{sb} proteins with locus tags: {len(pb)}\n")
    lg.write(f"reciprocal best hits: {n_rbh}  (ambiguous protein->locus: {n_ambig})\n")
    lg.write(f"ortholog rows written: {len(rows)}\n")
    lg.write(f"filters: pident>={min_id} qcoverage>={min_cov}\n")