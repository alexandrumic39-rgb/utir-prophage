#!/usr/bin/env python3
"""The ad-hoc analyses of 31 October, collected so they are reproducible.

Each of these was originally run as an inline shell heredoc. The results are
recorded in the project log; this file exists so they can be regenerated.

Run from ~/utir. Usage:
    python3 utir_adhoc.py            # all
    python3 utir_adhoc.py c e        # selected, by letter
"""

import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

GEN = "results/genomad_v113"
SEQ = "results/prophage_seqs_v113"
PHK = "results/pharokka_v113"
TAX = "results/gtdbtk_all/gtdbtk.bac120.summary.tsv"
CLS = "results/element_classification_v113.tsv"
CLUST = "results/isolate_clusters.tsv"
PANI = "results/prophage_ani_v113.tsv"
CHECKV = "results/checkv/quality_summary.tsv"
PRO = "results/prophages_per_isolate_v113.tsv"

VOTU_ANI, VOTU_AF, LEN_RATIO = 95.0, 85.0, 0.70


# ---------------------------------------------------------------- helpers

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


def verdicts():
    v = {}
    for r in csv.DictReader(open(CLS), delimiter="\t"):
        v[r["element"].replace("|", "__")] = r
    return v


def seq_lengths():
    return {os.path.basename(f)[:-4]:
            sum(len(l.strip()) for l in open(f) if not l.startswith(">"))
            for f in glob.glob(f"{SEQ}/*.fna")}


def ani_rows(path):
    """Yield (name_a, name_b, ani, af_ref, af_query) from a skani sparse file."""
    if not os.path.exists(path):
        return
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or f[0].startswith("Ref"):
            continue
        a, b = os.path.basename(f[0]), os.path.basename(f[1])
        a = a[:-4] if a.endswith(".fna") else a
        b = b[:-4] if b.endswith(".fna") else b
        try:
            yield a, b, float(f[2]), float(f[3]), float(f[4])
        except ValueError:
            continue


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs) ** 0.5
    vy = sum((y - my) ** 2 for y in ys) ** 0.5
    return cov / (vx * vy) if vx and vy else float("nan")


# ---------------------------------------------------------------- analyses

def a_sampling_structure():
    """Participant and timepoint structure, derived from filenames."""
    print("\n=== A. sampling structure from filenames ===")
    byp = defaultdict(set)
    for f in sorted(glob.glob("data/*.fna")):
        b = os.path.basename(f)[:-4]
        m = re.match(r"utir_(\d+)_([a-z0-9]+)_\d+$", b)
        if m:
            byp[m.group(1)].add(m.group(2))
        else:
            byp[b] = set()
    longitudinal = {p: t for p, t in byp.items() if len(t) > 1}
    print(f"  participants parsed            : {len(byp)}")
    print(f"  sampled at >1 timepoint        : {len(longitudinal)}")
    for p, t in sorted(longitudinal.items(), key=lambda kv: -len(kv[1])):
        print(f"    {p}: {len(t):>2}  {' '.join(sorted(t))}")
    single = sorted(p for p, t in byp.items() if len(t) == 1)
    print(f"  single-timepoint participants  : {len(single)}")
    print(f"    {' '.join(single)}")
    unparsed = sorted(p for p, t in byp.items() if not t)
    if unparsed:
        print(f"  did not match the naming scheme: {' '.join(unparsed)}")


def b_edge_by_verdict():
    """Contig-edge flag cross-tabulated against element verdict."""
    print("\n=== B. contig-edge flag by element verdict ===")
    edge = {}
    for f in glob.glob(f"{GEN}/*_find_proviruses/*_provirus.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            edge[r["seq_name"].replace("|", "__")] = r.get("in_seq_edge")
    tally = defaultdict(lambda: [0, 0])
    for el, r in verdicts().items():
        t = tally[r["verdict"]]
        t[0] += 1
        if edge.get(el) == "True":
            t[1] += 1
    print(f"  {'verdict':<16}{'n':>6}{'at edge':>10}")
    for v, (n, k) in sorted(tally.items()):
        print(f"  {v:<16}{n:>6}{k:>7} ({100 * k / n:.0f}%)")
    ends = sum(2 * sum(1 for l in open(f) if l.startswith(">"))
               for f in glob.glob("data/*.fna"))
    total = len(edge)
    at = sum(1 for v in edge.values() if v == "True")
    print(f"\n  contig ends in the dataset : {ends}")
    print(f"  provirus candidates        : {total}")
    print(f"  flagged at a contig edge   : {at} ({100 * at / total:.0f}% of calls)")


def c_votu_chaining():
    """Length heterogeneity within vOTUs, under three clustering rules."""
    print("\n=== C. vOTU clustering: chaining diagnostic ===")
    L = seq_lengths()
    if not L:
        print("  no prophage sequences found")
        return

    def cluster(rule):
        parent = {p: p for p in L}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b, ani, afr, afq in ani_rows(PANI):
            if a not in parent or b not in parent:
                continue
            if ani < VOTU_ANI:
                continue
            if rule == "min" and min(afr, afq) < VOTU_AF:
                continue
            if rule in ("max", "guard") and max(afr, afq) < VOTU_AF:
                continue
            if rule == "guard":
                lo, hi = min(L[a], L[b]), max(L[a], L[b])
                if hi and lo / hi < LEN_RATIO:
                    continue
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        g = defaultdict(list)
        for p in L:
            g[find(p)].append(p)
        return list(g.values())

    for rule, label in (("min", "min(AF) >= 85"),
                        ("max", "max(AF) >= 85"),
                        ("guard", f"max(AF) >= 85 and ratio >= {LEN_RATIO}")):
        groups = cluster(rule)
        singles = sum(1 for m in groups if len(m) == 1)
        bad = [(len(m), min(L[x] for x in m), max(L[x] for x in m))
               for m in groups
               if len(m) > 1 and min(L[x] for x in m) / max(L[x] for x in m) < 0.5]
        print(f"\n  {label:<40} vOTUs={len(groups):>4}  singletons={singles:>4}"
              f"  spanning >2x in length={len(bad):>3}")
        for n, lo, hi in sorted(bad, key=lambda t: t[1] / t[2])[:5]:
            print(f"      {n:>3} members  {lo:>8,}-{hi:>8,} bp  ratio {lo / hi:.2f}")


def d_whole_contig_twins():
    """Whole viral contigs that duplicate an integrated provirus."""
    print("\n=== D. whole viral contigs with a chromosomal twin ===")
    whole = defaultdict(list)
    for f in sorted(glob.glob(f"{GEN}/*_summary/*_virus_summary.tsv")):
        iso = os.path.basename(os.path.dirname(f))[: -len("_summary")]
        for r in csv.DictReader(open(f), delimiter="\t"):
            if "|provirus_" not in r["seq_name"]:
                whole[iso].append(r["seq_name"])
    twin = defaultdict(list)
    for a, b, ani, afr, afq in ani_rows(PANI):
        if iso_of(a) != iso_of(b) or ani < 95 or max(afr, afq) < 90:
            continue
        wa, wb = "__provirus_" not in a, "__provirus_" not in b
        if wa == wb:
            continue
        w, pv = (a, b) if wa else (b, a)
        twin[w].append((pv, ani))
    sp, wt, nt = species(), 0, 0
    for iso in sorted(whole):
        p = f"data/{iso}.fna"
        L = fasta_lengths(p) if os.path.exists(p) else {}
        for c in whole[iso]:
            t = twin.get(c, [])
            wt, nt = (wt + 1, nt) if t else (wt, nt + 1)
            tag = "twinned " if t else "isolated"
            print(f"  {tag} {c:<30}{L.get(c, 0):>9,} bp  {sp.get(iso, '?')[:28]}")
            for pv, ani in t[:2]:
                print(f"              twin {pv.split('__')[-1]}  ANI {ani:.1f}")
    print(f"\n  twinned (induction candidates)        : {wt}")
    print(f"  isolated (episome or assembly break)  : {nt}")


def e_fragmentation():
    """Contig count against prophage count, and genome size against load."""
    print("\n=== E. fragmentation and genome size versus prophage count ===")
    contigs, gb = {}, {}
    for f in glob.glob("data/*.fna"):
        iso = os.path.basename(f)[:-4]
        L = fasta_lengths(f)
        contigs[iso], gb[iso] = len(L), sum(L.values())
    pc = {}
    for line in open(PRO):
        a = line.rstrip("\n").split("\t")
        if len(a) == 2:
            pc[a[0]] = int(a[1])
    sp = species()
    iso = [i for i in pc if i in contigs]
    print(f"  {'prophages':>9}{'contigs':>9}{'Mb':>7}  isolate / species")
    for i in sorted(iso, key=lambda x: -pc[x])[:12]:
        print(f"  {pc[i]:>9}{contigs[i]:>9}{gb[i] / 1e6:>7.2f}  {i}  "
              f"{sp.get(i, '?')[:26]}")
    ys = [pc[i] for i in iso]
    print(f"\n  r (contig count vs prophage count) : "
          f"{pearson([contigs[i] for i in iso], ys):.3f}")
    print(f"  r (genome size  vs prophage count) : "
          f"{pearson([gb[i] / 1e6 for i in iso], ys):.3f}")
    dist = Counter(contigs[i] for i in iso)
    print(f"  contig count distribution          : "
          + "  ".join(f"{k}:{v}" for k, v in sorted(dist.items())))


def f_mobileog_tally():
    """mobileOG-db categories for the non-phage elements."""
    print("\n=== F. mobileOG categories for the non-phage elements ===")
    p = "results/mobileog_nonphage.tsv"
    if not os.path.exists(p):
        print(f"  {p} not found; run the DIAMOND search first")
        return
    major, gene, per = Counter(), Counter(), defaultdict(Counter)
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        s = f[1].split("|")
        cat = s[3] if len(s) > 3 else "?"
        major[cat] += 1
        if len(s) > 1:
            gene[s[1]] += 1
        per[re.sub(r"_CDS_\d+$", "", f[0])][cat] += 1
    tot = sum(major.values())
    print(f"  hits: {tot}   elements with a hit: {len(per)}")
    print(f"\n  {'major category':<36}{'hits':>6}{'share':>8}")
    for k, v in major.most_common():
        print(f"  {k:<36}{v:>6}{100 * v / tot:>7.1f}%")
    print(f"\n  most frequent gene names:")
    for k, v in gene.most_common(15):
        print(f"    {v:>4}  {k}")
    print("\n  Conjugation machinery is what distinguishes an ICE from a satellite:"
          "\n  a low 'transfer' share argues for satellites.")


def g_checkv_crosstab():
    """CheckV quality and contamination against element verdict."""
    print("\n=== G. CheckV quality and contamination by verdict ===")
    if not os.path.exists(CHECKV):
        print(f"  {CHECKV} not found; run checkv end_to_end first")
        return
    cv = {r["contig_id"]: r
          for r in csv.DictReader(open(CHECKV), delimiter="\t")}
    vd = {k: r["verdict"] for k, r in verdicts().items()}
    tab, cont = defaultdict(Counter), defaultdict(list)
    for cid, r in cv.items():
        v = vd.get(cid, "not classified")
        tab[v][r["checkv_quality"]] += 1
        try:
            cont[v].append(float(r["contamination"]))
        except (TypeError, ValueError):
            pass
    qs = ["Complete", "High-quality", "Medium-quality", "Low-quality",
          "Not-determined"]
    print(f"  {'verdict':<16}" + "".join(f"{q[:9]:>11}" for q in qs)
          + f"{'mean cont%':>11}")
    for v in sorted(tab):
        row = "".join(f"{tab[v].get(q, 0):>11}" for q in qs)
        c = cont[v]
        m = f"{sum(c) / len(c):>11.1f}" if c else f"{'-':>11}"
        print(f"  {v:<16}{row}{m}")
    allc = sorted(x for c in cont.values() for x in c)
    n = len(allc)
    if n:
        print(f"\n  contamination across {n} elements: median {allc[n // 2]:.1f}%,"
              f" mean {sum(allc) / n:.1f}%,"
              f" >25%: {sum(1 for x in allc if x > 25)},"
              f" >50%: {sum(1 for x in allc if x > 50)}")
    overall = Counter(r["checkv_quality"] for r in cv.values())
    print(f"\n  overall quality: "
          + "  ".join(f"{k} {v}" for k, v in overall.most_common()))


def h_strain_splits():
    """Same species and timepoint, different strain cluster."""
    print("\n=== H. polyclonal colonisation within a single sample ===")
    rows = list(csv.DictReader(open(CLUST), delimiter="\t"))
    sp = species()
    ani = {frozenset((a, b)): (x, y, z) for a, b, x, y, z in ani_rows(GANI_PATH)}
    grp = defaultdict(list)
    for r in rows:
        if r["participant_id"] == "NA":
            continue
        grp[(r["participant_id"], r["timepoint"], r["species_cluster"])].append(r)
    found = 0
    for k, g in sorted(grp.items()):
        if len(g) < 2 or len({r["strain_cluster"] for r in g}) < 2:
            continue
        found += 1
        print(f"\n  participant {k[0]}  {k[1]}  species cluster {k[2]}"
              f"  ({sp.get(g[0]['isolate_id'], '?')})")
        for r in g:
            print(f"    {r['isolate_id']:<24} strain {r['strain_cluster']}")
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                a, b = g[i]["isolate_id"], g[j]["isolate_id"]
                v = ani.get(frozenset((a, b)))
                print(f"    {a} vs {b}: "
                      + (f"ANI {v[0]}  AF {v[1]}/{v[2]}" if v
                         else "no ANI edge above threshold"))
    print(f"\n  cases found: {found}")


GANI_PATH = "results/ani_edges.tsv"

ANALYSES = {
    "a": a_sampling_structure,
    "b": b_edge_by_verdict,
    "c": c_votu_chaining,
    "d": d_whole_contig_twins,
    "e": e_fragmentation,
    "f": f_mobileog_tally,
    "g": g_checkv_crosstab,
    "h": h_strain_splits,
}

if __name__ == "__main__":
    wanted = [x.lower() for x in sys.argv[1:]] or sorted(ANALYSES)
    for key in wanted:
        fn = ANALYSES.get(key)
        if fn is None:
            print(f"\n!! no analysis '{key}'; available: {' '.join(sorted(ANALYSES))}")
            continue
        try:
            fn()
        except Exception as exc:
            print(f"\n!! {fn.__name__} failed: {type(exc).__name__}: {exc}")
