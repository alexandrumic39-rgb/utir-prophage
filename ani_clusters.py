#!/usr/bin/env python3
"""Cluster UTIr isolates by ANI and report recurrence within participants."""
import os
from collections import defaultdict

EDGES, DATA, OUT = "results/ani_edges.tsv", "data", "results/isolate_clusters.tsv"
SPECIES_ANI, SPECIES_AF = 95.0, 50.0
STRAIN_ANI, STRAIN_AF = 99.5, 80.0

def iso_name(path):
    b = os.path.basename(path)
    return b[:-4] if b.endswith(".fna") else b

def parse_id(iid):
    p = iid.split("_")
    if len(p) == 4 and p[1].isdigit():
        return p[1], p[2], p[3]
    return None, None, None

def components(nodes, edges):
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    groups = defaultdict(list)
    for n in nodes: groups[find(n)].append(n)
    labels = {}
    for i, members in enumerate(sorted(groups.values(), key=lambda m: (-len(m), sorted(m)[0])), 1):
        for m in members: labels[m] = i
    return labels

def main():
    isolates = sorted(f[:-4] for f in os.listdir(DATA) if f.endswith(".fna"))
    sp_edges, st_edges = [], []
    with open(EDGES) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5: continue
            try: ani, afr, afq = float(f[2]), float(f[3]), float(f[4])
            except ValueError: continue
            a, b, af = iso_name(f[0]), iso_name(f[1]), min(afr, afq)
            if ani >= SPECIES_ANI and af >= SPECIES_AF: sp_edges.append((a, b))
            if ani >= STRAIN_ANI and af >= STRAIN_AF: st_edges.append((a, b))
    sp, st = components(isolates, sp_edges), components(isolates, st_edges)
    with open(OUT, "w") as out:
        out.write("isolate_id\tparticipant_id\ttimepoint\tisolate_no\tspecies_cluster\tstrain_cluster\n")
        for i in isolates:
            p, tp, no = parse_id(i)
            out.write(f"{i}\t{p or 'NA'}\t{tp or 'NA'}\t{no or 'NA'}\t{sp[i]}\t{st[i]}\n")
    print(f"isolates                 : {len(isolates)}")
    print(f"species clusters         : {len(set(sp.values()))}   (ANI>={SPECIES_ANI}, AF>={SPECIES_AF}%)")
    print(f"strain clusters          : {len(set(st.values()))}   (ANI>={STRAIN_ANI}, AF>={STRAIN_AF}%)")
    sizes = defaultdict(int)
    for v in sp.values(): sizes[v] += 1
    big = sorted(sizes.values(), reverse=True)
    print(f"species cluster sizes    : {big}")
    print(f"singleton species        : {sum(1 for s in big if s == 1)}\n")
    for label, assign in (("SPECIES", sp), ("STRAIN", st)):
        byp = defaultdict(lambda: defaultdict(set))
        for i in isolates:
            p, tp, _ = parse_id(i)
            if p: byp[p][assign[i]].add(tp)
        print(f"=== same {label} at >=2 timepoints in one participant ===")
        found = False
        for p in sorted(byp):
            for c, tps in sorted(byp[p].items()):
                if len(tps) >= 2:
                    found = True
                    print(f"  participant {p}  cluster {c:>3}  timepoints: {' '.join(sorted(tps))}")
        if not found: print("  none")
        print()
    cl_parts = defaultdict(set)
    for i in isolates:
        p, _, _ = parse_id(i)
        if p: cl_parts[sp[i]].add(p)
    shared = {c: ps for c, ps in cl_parts.items() if len(ps) > 1}
    print(f"=== species clusters present in >1 participant: {len(shared)} ===")
    for c, ps in sorted(shared.items(), key=lambda kv: -len(kv[1]))[:12]:
        print(f"  cluster {c:>3}: {len(ps)} participants  ({' '.join(sorted(ps))})")
    print(f"\nwritten: {OUT}")

if __name__ == "__main__":
    main()
