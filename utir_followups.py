#!/usr/bin/env python3
"""Follow-up analyses on the UTIr prophage catalogue, from existing data."""

import csv
import glob
import os
import re
from collections import defaultdict

GEN = "results/genomad_v113"
SEQ = "results/prophage_seqs_v113"
PHK = "results/pharokka_v113"
TAX = "results/gtdbtk_all/gtdbtk.bac120.summary.tsv"
PRO = "results/prophages_per_isolate_v113.tsv"
CLS = "results/element_classification_v113.tsv"


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


def prophage_counts():
    d = {}
    for line in open(PRO):
        f = line.rstrip("\n").split("\t")
        if len(f) == 2:
            d[f[0]] = int(f[1])
    return d


def genome_sizes():
    out = {}
    for f in glob.glob("data/*.fna"):
        L = fasta_lengths(f)
        out[os.path.basename(f)[:-4]] = (sum(L.values()), len(L))
    return out


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs) ** 0.5
    vy = sum((y - my) ** 2 for y in ys) ** 0.5
    return cov / (vx * vy) if vx and vy else float("nan")


def a_load_per_mb():
    print("\n=== 1. prophage load per megabase, by species (n>=2) ===")
    sp, pc, gs = species(), prophage_counts(), genome_sizes()
    agg = defaultdict(lambda: [0, 0, 0.0])
    for iso, s in sp.items():
        if iso not in pc or iso not in gs:
            continue
        a = agg[s]
        a[0] += 1
        a[1] += pc[iso]
        a[2] += gs[iso][0] / 1e6
    print(f"{'species':<38}{'n':>3}{'Mb':>7}{'mean':>7}{'perMb':>7}")
    rows = [(k, v) for k, v in agg.items() if v[0] >= 2]
    for k, (n, t, mb) in sorted(rows, key=lambda kv: -kv[1][1] / kv[1][2]):
        print(f"{k:<38}{n:>3}{mb / n:>7.2f}{t / n:>7.2f}{t / mb:>7.2f}")
    iso = [i for i in pc if i in gs]
    xs = [gs[i][0] / 1e6 for i in iso]
    ys = [pc[i] for i in iso]
    print(f"\n  r (genome size vs prophage count), n={len(xs)}: {pearson(xs, ys):.3f}")
    per = [pc[i] / (gs[i][0] / 1e6) for i in iso]
    print(f"  prophages per Mb: mean {sum(per) / len(per):.2f}, "
          f"min {min(per):.2f}, max {max(per):.2f}")


def b_sizes():
    print("\n=== 2. prophage size distribution (>=5 kb) ===")
    L = sorted(sum(len(l.strip()) for l in open(f) if not l.startswith(">"))
               for f in glob.glob(f"{SEQ}/*.fna"))
    n = len(L)
    if not n:
        print("  no sequences found")
        return

    def q(p):
        return L[int(p * (n - 1))]

    print(f"  n={n}  min={L[0]:,}  Q1={q(.25):,}  median={q(.5):,}  "
          f"Q3={q(.75):,}  max={L[-1]:,}  mean={sum(L) // n:,}")
    for lo, hi in [(5, 10), (10, 20), (20, 30), (30, 40), (40, 50),
                   (50, 75), (75, 100), (100, 1000)]:
        c = sum(1 for x in L if lo * 1000 <= x < hi * 1000)
        print(f"  {lo:>4}-{hi:<4} kb {c:>4}  {'#' * (c // 3)}")


def c_terl():
    print("\n=== 3. elements carrying a terminase large subunit ===")
    p = f"{PHK}/terL.faa"
    if not os.path.exists(p):
        print("  terL.faa not found")
        return
    have = set()
    for line in open(p):
        if line.startswith(">"):
            have.add(re.sub(r"_CDS_\d+$", "", line[1:].split()[0]))
    tally = defaultdict(lambda: [0, 0])
    rows = list(csv.DictReader(open(CLS), delimiter="\t"))
    for r in rows:
        key = r["element"].replace("|", "__")
        t = tally[r["verdict"]]
        t[0] += 1
        if key in have:
            t[1] += 1
    print(f"  {'verdict':<16}{'n':>5}{'with terL':>12}")
    for v, (n, k) in sorted(tally.items()):
        print(f"  {v:<16}{n:>5}{k:>8} ({100 * k / n:.0f}%)")
    tot = sum(t[1] for t in tally.values())
    print(f"\n  total with terL: {tot} of {len(rows)} elements")


def d_cargo():
    for label, fn, hit, idc in (
            ("virulence (VFDB)", "top_hits_vfdb.tsv", "vfdb_hit", "vfdb_seqIdentity"),
            ("resistance (CARD)", "top_hits_card.tsv", "card_hit", "card_seqIdentity")):
        print(f"\n=== 4. cargo genes: {label} ===")
        p = f"{PHK}/{fn}"
        if not os.path.exists(p):
            print("  missing")
            continue
        for row in csv.DictReader(open(p), delimiter="\t"):
            c = row.get("contig", "")
            print(f"  {c[:44]:<44} {row.get(idc, ''):>6}  {row.get(hit, '')[:95]}")


def e_whole_contigs():
    print("\n=== 5. whole viral contigs in their genome context ===")
    calls = defaultdict(list)
    for f in sorted(glob.glob(f"{GEN}/*_summary/*_virus_summary.tsv")):
        iso = os.path.basename(os.path.dirname(f))[: -len("_summary")]
        for row in csv.DictReader(open(f), delimiter="\t"):
            if "|provirus_" not in row["seq_name"]:
                calls[iso].append(row["seq_name"])
    for iso in sorted(calls):
        path = f"data/{iso}.fna"
        if not os.path.exists(path):
            continue
        L = fasta_lengths(path)
        tot = sum(L.values())
        order = sorted(L.items(), key=lambda kv: -kv[1])
        rank = {k: i + 1 for i, (k, _) in enumerate(order)}
        print(f"\n  {iso}: {tot / 1e6:.2f} Mb in {len(L)} contigs")
        for c in calls[iso]:
            n = L.get(c, 0)
            print(f"    {c:<32}{n:>10,} bp  rank {rank.get(c, '?')}/{len(L)}"
                  f"  {100 * n / tot:>5.1f}% of genome")


def f_edges():
    print("\n=== 6. contig ends available vs occupied by a prophage ===")
    ends = 0
    for f in glob.glob("data/*.fna"):
        ends += 2 * sum(1 for l in open(f) if l.startswith(">"))
    total = at_edge = 0
    for f in glob.glob(f"{GEN}/*_find_proviruses/*_provirus.tsv"):
        for row in csv.DictReader(open(f), delimiter="\t"):
            total += 1
            if row.get("in_seq_edge") == "True":
                at_edge += 1
    print(f"  contig ends in the dataset       : {ends}")
    print(f"  provirus calls                   : {total}")
    print(f"  calls flagged at a contig edge   : {at_edge}")
    if ends:
        print(f"  share of contig ends occupied    : {100 * at_edge / ends:.1f}%")


if __name__ == "__main__":
    for fn in (a_load_per_mb, b_sizes, c_terl, d_cargo, e_whole_contigs, f_edges):
        try:
            fn()
        except Exception as exc:
            print(f"\n!! {fn.__name__} failed: {type(exc).__name__}: {exc}")
