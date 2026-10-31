#!/usr/bin/env python3
"""Second round of follow-ups: remaining analyses and open discrepancies."""

import csv
import glob
import os
import re
from collections import Counter, defaultdict

GEN = "results/genomad_v113"
SEQ = "results/prophage_seqs_v113"
PHK = "results/pharokka_v113"
TAX = "results/gtdbtk_all/gtdbtk.bac120.summary.tsv"
CLS = "results/element_classification_v113.tsv"
PANI = "results/prophage_ani_v113.tsv"
GANI = "results/ani_edges.tsv"


def fasta_lengths(path):
    d, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            d[name] = 0
        elif name is not None:
            d[name] += len(line.strip())
    return d


def species():
    s = {}
    for row in csv.DictReader(open(TAX), delimiter="\t"):
        t = row["classification"].split(";")
        s[row["user_genome"]] = t[-1].replace("s__", "").strip() or "unclassified"
    return s


def iso_of(pid):
    return re.sub(r"_\d+$", "", re.sub(r"__provirus_.*$", "", pid))


def a_phage_taxonomy():
    print("\n=== A. prophage taxonomy (INPHARED mash hits) ===")
    p = f"{PHK}/pharokka_top_hits_mash_inphared.tsv"
    if not os.path.exists(p):
        print("  file not found")
        return
    with open(p) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = rdr.fieldnames or []
        rows = list(rdr)
    print(f"  rows: {len(rows)}")
    print(f"  columns: {', '.join(cols)}")
    shown = 0
    for key in ("Genus", "genus", "Family", "family", "Host", "host",
                "Description", "description", "Classification"):
        if key in cols:
            t = Counter((r.get(key) or "unassigned").strip() for r in rows)
            print(f"\n  tally of '{key}' ({len(t)} distinct):")
            for k, v in t.most_common(25):
                print(f"    {v:>4}  {k[:78]}")
            shown += 1
            if shown == 2:
                break
    if not shown and rows:
        print("  no recognised taxonomic column; first row:")
        for k, v in list(rows[0].items())[:14]:
            print(f"    {k}: {str(v)[:88]}")


def b_whole_contig_split():
    print("\n=== B. whole viral contigs: with vs without a chromosomal twin ===")
    whole = defaultdict(list)
    for f in sorted(glob.glob(f"{GEN}/*_summary/*_virus_summary.tsv")):
        iso = os.path.basename(os.path.dirname(f))[: -len("_summary")]
        for row in csv.DictReader(open(f), delimiter="\t"):
            if "|provirus_" not in row["seq_name"]:
                whole[iso].append(row["seq_name"])
    twin = defaultdict(list)
    if os.path.exists(PANI):
        for line in open(PANI):
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0].startswith("Ref"):
                continue
            a = os.path.basename(f[0])[:-4]
            b = os.path.basename(f[1])[:-4]
            try:
                ani, afr, afq = float(f[2]), float(f[3]), float(f[4])
            except ValueError:
                continue
            if iso_of(a) != iso_of(b) or ani < 95 or max(afr, afq) < 90:
                continue
            wa, wb = "__provirus_" not in a, "__provirus_" not in b
            if wa == wb:
                continue
            w, pv = (a, b) if wa else (b, a)
            twin[w].append((pv, ani))
    sp, nt, wt = species(), 0, 0
    for iso in sorted(whole):
        path = f"data/{iso}.fna"
        L = fasta_lengths(path) if os.path.exists(path) else {}
        for c in whole[iso]:
            t = twin.get(c, [])
            if t:
                wt += 1
                tag = "INDUCTION?"
            else:
                nt += 1
                tag = "episome?  "
            print(f"  {tag} {c:<30}{L.get(c, 0):>9,} bp  {sp.get(iso, '?')[:30]}")
            for pv, ani in t[:2]:
                print(f"               twin: {pv.split('__')[-1]}  ANI {ani:.1f}")
    print(f"\n  with a chromosomal twin (induction candidates)     : {wt}")
    print(f"  without a twin (plasmid-like prophage candidates)  : {nt}")


def c_large_elements():
    print("\n=== C. elements of 75 kb and above ===")
    big = []
    for f in glob.glob(f"{SEQ}/*.fna"):
        n = sum(len(l.strip()) for l in open(f) if not l.startswith(">"))
        if n >= 75000:
            big.append((n, os.path.basename(f)[:-4]))
    cls = {}
    for r in csv.DictReader(open(CLS), delimiter="\t"):
        cls[r["element"].replace("|", "__")] = r
    sp = species()
    for n, pid in sorted(big, reverse=True):
        r = cls.get(pid, {})
        print(f"\n  {pid}")
        print(f"    {n:,} bp   {sp.get(iso_of(pid), '?')}")
        print(f"    verdict={r.get('verdict', '?')}  genes={r.get('n_genes', '?')}"
              f"  hallmarks={r.get('hallmarks', '?')}"
              f"  structural={r.get('structural', '?')}"
              f"  integrase={r.get('integrase', '?')}")


def d_reconcile():
    print("\n=== D. reconciling provirus candidates against final calls ===")
    cand = set()
    for f in glob.glob(f"{GEN}/*_find_proviruses/*_provirus.tsv"):
        for row in csv.DictReader(open(f), delimiter="\t"):
            cand.add(row["seq_name"])
    final = set()
    for f in glob.glob(f"{GEN}/*_summary/*_virus_summary.tsv"):
        for row in csv.DictReader(open(f), delimiter="\t"):
            final.add(row["seq_name"])
    fin_pv = {s for s in final if "|provirus_" in s}
    print(f"  provirus candidates from find-proviruses : {len(cand)}")
    print(f"  final calls in the summaries             : {len(final)}")
    print(f"    integrated proviruses                  : {len(fin_pv)}")
    print(f"    whole viral contigs                    : {len(final) - len(fin_pv)}")
    dropped = sorted(cand - fin_pv)
    print(f"  candidates dropped by the classifier     : {len(dropped)}")
    for s in dropped[:25]:
        print(f"    {s}")
    extra = sorted(fin_pv - cand)
    if extra:
        print(f"  final proviruses absent from candidates  : {len(extra)}")


def genome_ani():
    d = {}
    if not os.path.exists(GANI):
        return d
    for line in open(GANI):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[0].startswith("Ref"):
            continue
        a, b = os.path.basename(f[0]), os.path.basename(f[1])
        a = a[:-4] if a.endswith(".fna") else a
        b = b[:-4] if b.endswith(".fna") else b
        d[frozenset((a, b))] = (f[2], f[3], f[4])
    return d


def e_strain_splits():
    print("\n=== E. same species and timepoint, different strain cluster ===")
    rows = list(csv.DictReader(open("results/isolate_clusters.tsv"), delimiter="\t"))
    sp, ani = species(), genome_ani()
    grp = defaultdict(list)
    for r in rows:
        if r["participant_id"] == "NA":
            continue
        grp[(r["participant_id"], r["timepoint"], r["species_cluster"])].append(r)
    found = 0
    for k, g in sorted(grp.items()):
        strains = {r["strain_cluster"] for r in g}
        if len(g) < 2 or len(strains) < 2:
            continue
        found += 1
        print(f"\n  participant {k[0]}  {k[1]}  species cluster {k[2]}"
              f"  ({sp.get(g[0]['isolate_id'], '?')})")
        for r in g:
            print(f"    {r['isolate_id']:<24} strain cluster {r['strain_cluster']}")
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                a, b = g[i]["isolate_id"], g[j]["isolate_id"]
                v = ani.get(frozenset((a, b)))
                if v:
                    print(f"    {a} vs {b}: ANI {v[0]}  AF {v[1]}/{v[2]}")
                else:
                    print(f"    {a} vs {b}: no ANI edge above threshold")
    if not found:
        print("  none")


def f_superseded():
    print("\n=== F. geNomad runs on disk, superseded and current ===")
    for d in sorted(glob.glob("results/genomad*")):
        if not os.path.isdir(d):
            continue
        summ = sorted(glob.glob(f"{d}/*_summary"))
        tot = 0
        for s in summ:
            iso = os.path.basename(s)[: -len("_summary")]
            p = f"{s}/{iso}_virus_summary.tsv"
            if os.path.exists(p):
                tot += max(0, sum(1 for _ in open(p)) - 1)
        print(f"  {d:<36} genomes={len(summ):>4}  calls={tot:>5}")


if __name__ == "__main__":
    for fn in (a_phage_taxonomy, b_whole_contig_split, c_large_elements,
               d_reconcile, e_strain_splits, f_superseded):
        try:
            fn()
        except Exception as exc:
            print(f"\n!! {fn.__name__} failed: {type(exc).__name__}: {exc}")
