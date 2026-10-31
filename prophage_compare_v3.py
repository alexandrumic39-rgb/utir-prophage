#!/usr/bin/env python3
"""Compare prophage content across timepoints within persisting strain lineages."""

import glob
import os
import subprocess
import sys
from collections import defaultdict

GENOMAD = "results/genomad_v113"
CLUSTERS = "results/isolate_clusters.tsv"
SEQDIR = "results/prophage_seqs_v113"
EDGES = "results/prophage_ani_v113.tsv"
MATRIX = "results/lineage_prophage_matrix_v113_guarded.tsv"

VOTU_ANI, VOTU_AF = 95.0, 85.0
LEN_RATIO = 0.70
MIN_LEN = 5000
THREADS = "16"


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def load_meta():
    meta = {}
    with open(CLUSTERS) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                meta[f[0]] = {"participant": f[1], "timepoint": f[2],
                              "isolate_no": f[3], "species": f[4], "strain": f[5]}
    return meta


def collect_prophages():
    os.makedirs(SEQDIR, exist_ok=True)
    for old in glob.glob(f"{SEQDIR}/*.fna"):
        os.remove(old)
    info = {}
    for summ in sorted(glob.glob(f"{GENOMAD}/*_summary")):
        iso = os.path.basename(summ)
        iso = iso[: -len("_summary")]
        fna = os.path.join(summ, f"{iso}_virus.fna")
        if not os.path.exists(fna):
            continue
        for header, seq in read_fasta(fna):
            if len(seq) < MIN_LEN:
                continue
            pid = header.split()[0].replace("|", "__").replace("/", "_")
            info[pid] = (iso, len(seq))
            with open(os.path.join(SEQDIR, pid + ".fna"), "w") as out:
                out.write(f">{pid}\n")
                for i in range(0, len(seq), 70):
                    out.write(seq[i:i + 70] + "\n")
    return info


def run_skani():
    files = sorted(glob.glob(f"{SEQDIR}/*.fna"))
    if len(files) < 2:
        open(EDGES, "w").write("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n")
        return
    subprocess.run(["skani", "triangle", "-t", THREADS, "--sparse",
                    "-o", EDGES] + files, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def votus(pids, slen):
    parent = {p: p for p in pids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    if os.path.exists(EDGES):
        with open(EDGES) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                try:
                    ani, afr, afq = float(f[2]), float(f[3]), float(f[4])
                except ValueError:
                    continue
                if ani < VOTU_ANI or max(afr, afq) < VOTU_AF:
                    continue
                a = os.path.basename(f[0])[:-4]
                b = os.path.basename(f[1])[:-4]
                if a in parent and b in parent:
                    if (min(slen[a], slen[b]) / max(slen[a], slen[b])
                            < LEN_RATIO):
                        continue
                    ra, rb = find(a), find(b)
                    if ra != rb:
                        parent[ra] = rb

    groups = defaultdict(list)
    for p in pids:
        groups[find(p)].append(p)
    labels = {}
    ordered = sorted(groups.values(), key=lambda m: (-len(m), sorted(m)[0]))
    for i, members in enumerate(ordered, 1):
        for m in members:
            labels[m] = i
    return labels


def main():
    meta = load_meta()
    info = collect_prophages()
    if not info:
        sys.exit("No prophage sequences found - has geNomad finished?")
    print(f"prophages >= {MIN_LEN} bp : {len(info)}")
    run_skani()
    votu = votus(sorted(info), {k: v[1] for k, v in info.items()})
    print(f"vOTUs (ANI>={VOTU_ANI}, AF>={VOTU_AF}%) : {len(set(votu.values()))}")

    by_iso = defaultdict(lambda: defaultdict(int))
    for pid, (iso, _) in info.items():
        by_iso[iso][votu[pid]] += 1

    analysed = sorted(os.path.basename(d)[: -len("_summary")]
                      for d in glob.glob(f"{GENOMAD}/*_summary"))
    lineages = defaultdict(list)
    for iso in analysed:
        m = meta.get(iso)
        if m and m["participant"] != "NA":
            lineages[(m["participant"], m["strain"])].append(iso)

    rows = []
    for (part, strain), isos in sorted(lineages.items()):
        tps = {meta[i]["timepoint"] for i in isos}
        if len(tps) < 2:
            continue
        isos = sorted(isos, key=lambda i: (meta[i]["timepoint"], meta[i]["isolate_no"]))
        present = sorted({v for i in isos for v in by_iso[i]})
        print()
        print(f"participant {part} | strain lineage {strain} | {len(isos)} isolates | {len(present)} vOTUs")
        print("  " + "vOTU".ljust(8) + "".join(
            f"{meta[i]['timepoint']}/{meta[i]['isolate_no']}".rjust(12) for i in isos))
        for v in present:
            print("  " + f"vOTU_{v}".ljust(8) + "".join(
                (str(by_iso[i][v]) if by_iso[i][v] else ".").rjust(12) for i in isos))
            for i in isos:
                rows.append((part, strain, meta[i]["timepoint"],
                             meta[i]["isolate_no"], i, v, by_iso[i][v]))
        for tp in sorted(tps):
            same = [i for i in isos if meta[i]["timepoint"] == tp]
            if len(same) < 2:
                continue
            sets = {tuple(sorted((k, c) for k, c in by_iso[i].items() if c)) for i in same}
            if len(sets) > 1:
                print(f"    note: colonies at {tp} disagree - pipeline noise, not change over time")

    with open(MATRIX, "w") as out:
        out.write("participant\tstrain_cluster\ttimepoint\tisolate_no\tisolate_id\tvotu\tcopies\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nwritten: {MATRIX}")


if __name__ == "__main__":
    main()
