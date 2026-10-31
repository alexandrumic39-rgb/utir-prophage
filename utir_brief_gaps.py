#!/usr/bin/env python3
"""Two analyses named in the project brief that need no new tools."""

import csv
import glob
import os
from collections import defaultdict

GEN = "results/genomad_v113"
TAX = "results/gtdbtk_all/gtdbtk.bac120.summary.tsv"
CLUST = "results/isolate_clusters.tsv"


def species():
    s = {}
    for row in csv.DictReader(open(TAX), delimiter="\t"):
        t = row["classification"].split(";")
        s[row["user_genome"]] = t[-1].replace("s__", "").strip() or "unclassified"
    return s


def genome_bp():
    out = {}
    for f in glob.glob("data/*.fna"):
        n = sum(len(l.strip()) for l in open(f) if not l.startswith(">"))
        out[os.path.basename(f)[:-4]] = n
    return out


def prophage_bp_and_count():
    bp, n = defaultdict(int), defaultdict(int)
    for d in sorted(glob.glob(f"{GEN}/*_summary")):
        iso = os.path.basename(d)[: -len("_summary")]
        p = os.path.join(d, f"{iso}_virus_summary.tsv")
        if not os.path.exists(p):
            continue
        n[iso] = 0
        with open(p) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for r in rows:
            n[iso] += 1
            for k in ("length", "seq_len", "Length"):
                if k in r:
                    try:
                        bp[iso] += int(r[k])
                    except (TypeError, ValueError):
                        pass
                    break
    return bp, n


def meta():
    m = {}
    for r in csv.DictReader(open(CLUST), delimiter="\t"):
        tp = r["timepoint"]
        if tp == "NA":
            kind = "unknown"
        elif tp.startswith("m"):
            kind = "UTI-free"
        elif tp.startswith("u"):
            kind = "UTI"
        else:
            kind = "unknown"
        m[r["isolate_id"]] = (r["participant_id"], tp, kind)
    return m


def a_uti_vs_free():
    print("\n=== 1. UTI versus UTI-free sampling ===")
    sp, gb, m = species(), genome_bp(), meta()
    pb, pc = prophage_bp_and_count()

    grp = defaultdict(list)
    for iso, (part, tp, kind) in m.items():
        if iso in pc and kind in ("UTI", "UTI-free"):
            grp[kind].append(iso)

    print(f"\n  {'group':<12}{'n':>4}{'carry':>7}{'mean':>7}{'perMb':>7}")
    for k in ("UTI-free", "UTI"):
        g = grp[k]
        if not g:
            continue
        tot = sum(pc[i] for i in g)
        mb = sum(gb[i] for i in g) / 1e6
        carry = sum(1 for i in g if pc[i])
        print(f"  {k:<12}{len(g):>4}{carry:>7}{tot / len(g):>7.2f}{tot / mb:>7.2f}")

    print("\n  Species composition differs between groups, so the raw means above"
          "\n  are not comparable. Composition:")
    comp = defaultdict(lambda: defaultdict(int))
    for k in ("UTI-free", "UTI"):
        for i in grp[k]:
            comp[sp.get(i, "?")][k] += 1
    print(f"  {'species':<38}{'free':>6}{'UTI':>6}")
    for s in sorted(comp, key=lambda x: -(comp[x]["UTI-free"] + comp[x]["UTI"])):
        c = comp[s]
        if c["UTI-free"] + c["UTI"] >= 3:
            print(f"  {s:<38}{c['UTI-free']:>6}{c['UTI']:>6}")

    print("\n  Within-species comparison, the only defensible view "
          "(species present in both groups):")
    print(f"  {'species':<34}{'free n':>7}{'free':>7}{'UTI n':>7}{'UTI':>7}")
    for s in sorted(comp):
        f = [i for i in grp["UTI-free"] if sp.get(i) == s]
        u = [i for i in grp["UTI"] if sp.get(i) == s]
        if f and u:
            mf = sum(pc[i] for i in f) / len(f)
            mu = sum(pc[i] for i in u) / len(u)
            print(f"  {s:<34}{len(f):>7}{mf:>7.2f}{len(u):>7}{mu:>7.2f}")

    print("\n  Per participant sampled in both states:")
    byp = defaultdict(lambda: defaultdict(list))
    for iso, (part, tp, kind) in m.items():
        if iso in pc and kind in ("UTI", "UTI-free"):
            byp[part][kind].append(iso)
    for p in sorted(byp):
        f, u = byp[p]["UTI-free"], byp[p]["UTI"]
        if f and u:
            mf = sum(pc[i] for i in f) / len(f)
            mu = sum(pc[i] for i in u) / len(u)
            print(f"    participant {p}: free n={len(f)} mean={mf:.2f}   "
                  f"UTI n={len(u)} mean={mu:.2f}")


def b_genomic_fraction():
    print("\n=== 2. prophage number and genomic fraction per isolate ===")
    sp, gb = species(), genome_bp()
    pb, pc = prophage_bp_and_count()
    rows = []
    for iso in pc:
        if iso in gb and gb[iso]:
            rows.append((100 * pb[iso] / gb[iso], pc[iso], pb[iso], gb[iso], iso))
    if not rows:
        print("  no data")
        return
    rows.sort(reverse=True)
    print(f"\n  highest prophage fraction:")
    print(f"  {'isolate':<22}{'calls':>6}{'prophage bp':>13}{'genome bp':>12}"
          f"{'%':>7}  species")
    for frac, n, p, g, iso in rows[:15]:
        print(f"  {iso:<22}{n:>6}{p:>13,}{g:>12,}{frac:>7.2f}  "
              f"{sp.get(iso, '?')[:26]}")
    fr = sorted(r[0] for r in rows)
    k = len(fr)
    print(f"\n  across {k} isolates: median {fr[k // 2]:.2f}%, "
          f"mean {sum(fr) / k:.2f}%, max {fr[-1]:.2f}%, min {fr[0]:.2f}%")
    agg = defaultdict(lambda: [0, 0.0])
    for frac, n, p, g, iso in rows:
        a = agg[sp.get(iso, "?")]
        a[0] += 1
        a[1] += frac
    print(f"\n  mean genomic fraction by species (n>=2):")
    print(f"  {'species':<40}{'n':>4}{'%':>8}")
    for s, (n, t) in sorted(agg.items(), key=lambda kv: -kv[1][1] / kv[1][0]):
        if n >= 2:
            print(f"  {s:<40}{n:>4}{t / n:>8.2f}")


if __name__ == "__main__":
    for fn in (a_uti_vs_free, b_genomic_fraction):
        try:
            fn()
        except Exception as exc:
            print(f"\n!! {fn.__name__} failed: {type(exc).__name__}: {exc}")
