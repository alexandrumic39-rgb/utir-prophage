#!/usr/bin/env python3
"""Decide which geNomad 'provirus' calls encode a real structural module."""
import glob, os, re, sys
from collections import defaultdict

GENOMAD = sys.argv[1] if len(sys.argv) > 1 else "results/genomad_long"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/element_classification.tsv"

STRUCTURAL = re.compile(r"capsid|portal|tail|terminase[, ]+large|major head|"
                        r"baseplate|tape measure|head-tail|neck|sheath|major tail", re.I)
INTEGRASE = re.compile(r"integrase|recombinase|transposase", re.I)
EXCISIONASE = re.compile(r"excisionase|\bxis\b|DUF1233", re.I)
LYSIS = re.compile(r"holin|endolysin|lysin|lysozyme|spanin", re.I)
GENE_SUFFIX = re.compile(r"_\d+$")

def main():
    elements = defaultdict(lambda: {"isolate": "", "genes": 0, "hallmarks": 0,
        "structural": 0, "integrase": 0, "excisionase": 0, "lysis": 0})
    files = sorted(glob.glob(f"{GENOMAD}/*_summary/*_virus_genes.tsv"))
    if not files:
        sys.exit(f"no *_virus_genes.tsv under {GENOMAD}/")
    for path in files:
        isolate = os.path.basename(os.path.dirname(path))[: -len("_summary")]
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                i_hall = header.index("virus_hallmark")
                i_desc = header.index("annotation_description")
            except ValueError:
                continue
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= i_desc:
                    continue
                rec = elements[GENE_SUFFIX.sub("", f[0])]
                rec["isolate"] = isolate
                rec["genes"] += 1
                try:
                    rec["hallmarks"] += int(f[i_hall])
                except ValueError:
                    pass
                d = f[i_desc]
                if d and d != "NA":
                    if STRUCTURAL.search(d): rec["structural"] += 1
                    if INTEGRASE.search(d): rec["integrase"] += 1
                    if EXCISIONASE.search(d): rec["excisionase"] += 1
                    if LYSIS.search(d): rec["lysis"] += 1

    def verdict(r):
        if r["structural"] >= 2 or (r["structural"] >= 1 and r["hallmarks"] >= 3):
            return "phage"
        if r["structural"] == 0 and (r["integrase"] or r["excisionase"]):
            return "satellite/ICE?"
        if r["structural"] == 0 and r["hallmarks"] <= 1:
            return "uncertain"
        return "phage-like"

    rows = [(r["isolate"], el, r["genes"], r["hallmarks"], r["structural"],
             r["integrase"], r["excisionase"], r["lysis"], verdict(r))
            for el, r in elements.items()]
    rows.sort(key=lambda x: (x[8], x[0], x[1]))
    with open(OUT, "w") as out:
        out.write("isolate\telement\tn_genes\thallmarks\tstructural\tintegrase\texcisionase\tlysis\tverdict\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    counts = defaultdict(int)
    for r in rows: counts[r[8]] += 1
    print(f"elements examined : {len(rows)}")
    for k in sorted(counts): print(f"  {k:<16} {counts[k]}")
    print("\nnon-phage candidates (no structural module):")
    print(f"  {'element':<52}{'genes':>6}{'hall':>5}{'str':>5}{'int':>5}{'xis':>5}")
    for r in rows:
        if r[8] in ("satellite/ICE?", "uncertain"):
            short = r[1].split("|")[-1] if "|" in r[1] else r[1]
            label = f"{r[0]} {short}"
            print(f"  {label[:52]:<52}{r[2]:>6}{r[3]:>5}{r[4]:>5}{r[5]:>5}{r[6]:>5}")
    print(f"\nwritten: {OUT}")

if __name__ == "__main__":
    main()
