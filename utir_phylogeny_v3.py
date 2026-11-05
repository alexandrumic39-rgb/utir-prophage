#!/usr/bin/env python3
"""
UTIr Prophage Phylogeny — v3
Fixes: species lookup (strips contig suffix), larger readable labels,
       better colour coding, cleaner layout.
Run from ~/utir:
    micromamba run -n pharokka18 python3 utir_phylogeny_v3.py
"""

import os, sys, csv, re, warnings
warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)
OUT_TREE = "results/utir_terl_phylogeny_v3.png"
OUT_DIST = "results/utir_terl_distances_v3.tsv"

REQUIRED = [
    "results/pharokka_v113/terL.faa",
    "results/gtdbtk_all/gtdbtk.bac120.summary.tsv",
    "results/pharokka_v113/top_hits_vfdb.tsv",
]
missing = [f for f in REQUIRED if not os.path.exists(f)]
if missing:
    print("ERROR: missing files:")
    for m in missing: print(f"  {m}")
    sys.exit(1)

print("Loading libraries...")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
print("  OK")

# =============================================================================
# HELPERS
# =============================================================================

def shorten_species(full):
    p = full.strip().split()
    return f"{p[0][0]}. {' '.join(p[1:])}" if len(p) >= 2 else full

def load_species():
    """GTDB-Tk → {isolate_id: 'E. coli'}
    Keys look like: utir_011_m0_01  (no contig suffix)
    """
    s = {}
    with open("results/gtdbtk_all/gtdbtk.bac120.summary.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genome = row["user_genome"]   # e.g. utir_011_m0_01
            parts  = row.get("classification","").split(";")
            sp     = parts[-1].replace("s__","").strip() or "unclassified"
            s[genome] = shorten_species(sp)
    print(f"  GTDB-Tk: {len(s)} isolates")
    # debug: show a few keys so we know the format
    sample = list(s.items())[:3]
    for k,v in sample:
        print(f"    example key: '{k}'  →  '{v}'")
    return s

def iso_from_pid(pid):
    """
    Convert a Pharokka terL sequence ID to its parent isolate ID.
    Examples:
      utir_011_m0_01_001__provirus_..._CDS_12  →  utir_011_m0_01
      utir_011_m0_01_002_CDS_5                 →  utir_011_m0_01
      utir_011_m0_01_001                       →  utir_011_m0_01
    Strategy: strip CDS suffix, strip provirus suffix, strip trailing _NNN contig number.
    """
    s = re.sub(r"_CDS_\d+$", "", pid)
    s = re.sub(r"(__|\|)provirus.*", "", s)
    # strip trailing contig number _001, _002 etc.
    s = re.sub(r"_\d{3}$", "", s)
    return s

def load_sit_isolates(species_map):
    """Return set of isolate IDs (matching species_map keys) carrying Sit genes."""
    sit = set()
    with open("results/pharokka_v113/top_hits_vfdb.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if "sit" in (row.get("vfdb_hit","") or "").lower():
                contig = row.get("contig","")
                iso    = iso_from_pid(contig)
                sit.add(iso)
    # verify against species map
    matched = {i for i in sit if i in species_map}
    unmatched = sit - matched
    print(f"  Sit isolates found: {len(sit)}  matched in GTDB: {len(matched)}")
    if unmatched:
        print(f"  Unmatched (will still flag): {sorted(unmatched)[:5]}")
    return sit   # keep all, matched or not

SHARED_ISOLATES = {
    "utir_025_m0_07","utir_029_m0_06",
    "utir_011_u1k2_02","utir_011_u1k2_04",
    "utir_006_m0_02","utir_010_u1k1_01","utir_010_u1k1_03",
}

# =============================================================================
# STEP 1 — LOAD SEQUENCES
# =============================================================================
print("\nStep 1: Loading terL sequences...")
species_map  = load_species()
sit_isolates = load_sit_isolates(species_map)

all_records = list(SeqIO.parse("results/pharokka_v113/terL.faa","fasta"))
print(f"  Total terL sequences: {len(all_records)}")

records_filtered = [r for r in all_records if len(r.seq) >= 200]
print(f"  After ≥200 aa filter: {len(records_filtered)}")

# keep longest per isolate
by_isolate = {}
for rec in records_filtered:
    iso = iso_from_pid(rec.id)
    if iso not in by_isolate or len(rec.seq) > len(by_isolate[iso][1].seq):
        by_isolate[iso] = (iso, rec)

records_use = [rec for iso, rec in by_isolate.values()]
print(f"  After one-per-isolate: {len(records_use)}")

# debug: show a few iso extractions
print("  ID extraction examples:")
for rec in records_use[:3]:
    iso = iso_from_pid(rec.id)
    sp  = species_map.get(iso, "NOT FOUND")
    print(f"    {rec.id[:45]}  →  iso='{iso}'  sp='{sp}'")

if len(records_use) > 60:
    records_use.sort(key=lambda r: -len(r.seq))
    records_use = records_use[:60]
    print(f"  Capped at 60 for readability")

# build annotations
annotations = {}   # short_id → (species, sit, shared, full_iso)
short_ids   = []

for rec in records_use:
    iso    = iso_from_pid(rec.id)
    sp     = species_map.get(iso, "unknown")
    sit    = iso in sit_isolates
    shared = iso in SHARED_ISOLATES
    # readable short label: drop 'utir_' prefix
    sid = iso.replace("utir_","")
    base, counter = sid, 1
    while sid in annotations:
        sid = f"{base}_{counter}"; counter += 1
    annotations[sid] = (sp, sit, shared, iso)
    short_ids.append(sid)

n_known = sum(1 for v in annotations.values() if v[0] != "unknown")
print(f"  Species resolved: {n_known}/{len(annotations)}")

# =============================================================================
# STEP 2 — ALIGNMENT
# =============================================================================
print("\nStep 2: Aligning sequences...")
max_len = max(len(r.seq) for r in records_use)
aligned = [SeqRecord(Seq(str(r.seq).ljust(max_len,"-")), id=sid, description="")
           for r, sid in zip(records_use, short_ids)]
alignment = MultipleSeqAlignment(aligned)
print(f"  {len(aligned)} sequences × {max_len} positions")

# =============================================================================
# STEP 3 — DISTANCES
# =============================================================================
print("\nStep 3: Calculating distances...")
calculator  = DistanceCalculator("blosum62")
dist_matrix = calculator.get_distance(alignment)

if   hasattr(dist_matrix, "ids"):   dm_ids = list(dist_matrix.ids)
elif hasattr(dist_matrix, "names"): dm_ids = list(dist_matrix.names)
else:                                dm_ids = short_ids
print(f"  {len(dm_ids)} × {len(dm_ids)} matrix")

# =============================================================================
# STEP 4 — NJ TREE
# =============================================================================
print("\nStep 4: Building NJ tree...")
tree = DistanceTreeConstructor().nj(dist_matrix)
tree.root_at_midpoint()
n_leaves = tree.count_terminals()
print(f"  {n_leaves} leaves")

# =============================================================================
# STEP 5 — PLOT
# =============================================================================
print("\nStep 5: Drawing tree...")

BG="#FFFFFF"; BORDER="#DEE2E6"; DARK="#111827"; MID="#374151"; SOFT="#9CA3AF"

SPECIES_COLS = {
    "E. coli":              "#DC2626",
    "E. faecalis":          "#D97706",
    "S. haemolyticus":      "#7C3AED",
    "S. vaginalis":         "#2563EB",
    "L. jensenii":          "#16A34A",
    "L. gasseri":           "#15803D",
    "E_B. lactis":          "#DB2777",
    "K. michiganensis":     "#0E7490",
    "P. mirabilis":         "#0284C7",
    "S. lugdunensis":       "#64748B",
    "S. epidermidis":       "#94A3B8",
    "A. mictus":            "#F59E0B",
    "W. anitrata":          "#CA8A04",
    "S. ureilyticus":       "#E11D48",
    "S. pettenkoferi":      "#F43F5E",
    "S. warneri":           "#10B981",
    "M. morganii":          "#1D4ED8",
    "P. rettgeri":          "#6D28D9",
    "C. freundii":          "#0891B2",
    "E. quasiroggenkampii": "#EA580C",
    "E. intestinihominis":  "#F97316",
    "unknown":              "#CBD5E1",
}

def sp_col(sp):
    if sp in SPECIES_COLS: return SPECIES_COLS[sp]
    for k,v in SPECIES_COLS.items():
        if sp.lower().startswith(k[:6].lower()): return v
    return "#CBD5E1"

# figure: wider annotation panel, taller rows
fig_h = max(14, n_leaves * 0.38)
fig, (ax_tree, ax_lab) = plt.subplots(
    1, 2, figsize=(20, fig_h), facecolor=BG,
    gridspec_kw={"width_ratios":[1.8, 1]})

# ── tree panel ────────────────────────────────────────────────────────────────
ax_tree.set_facecolor(BG)
Phylo.draw(tree, axes=ax_tree, do_show=False,
           label_func=lambda x:"", branch_labels=lambda x:"")
ax_tree.spines["top"].set_visible(False)
ax_tree.spines["right"].set_visible(False)
ax_tree.spines["left"].set_visible(False)
ax_tree.spines["bottom"].set_color(BORDER)
ax_tree.set_xlabel("Evolutionary distance (aa substitutions per site)",
                   color=MID, fontsize=10)
ax_tree.set_title(
    "Terminase Large Subunit (terL) Phylogeny — UTIr Prophage Catalogue\n"
    "Neighbour-Joining  ·  BLOSUM62  ·  midpoint rooted",
    color=DARK, fontsize=12, fontweight="bold", loc="left", pad=12)
ax_tree.set_yticks([])
ax_tree.tick_params(colors=MID, labelsize=9)

# ── annotation panel ──────────────────────────────────────────────────────────
ax_lab.set_facecolor(BG); ax_lab.axis("off")
term_order = [t.name for t in tree.get_terminals()]
n = len(term_order)
ax_lab.set_xlim(0, 1); ax_lab.set_ylim(0, n + 1)

# column headers
for x, lbl in [(0.01,"Isolate"), (0.42,"Species"), (0.80,"Flags")]:
    ax_lab.text(x, n+0.75, lbl, color=DARK, fontsize=9,
                fontweight="bold", va="center")
ax_lab.axhline(n+0.4, color=BORDER, linewidth=1)

fsz = max(6.5, min(9.5, 180/n))

for i, name in enumerate(reversed(term_order)):
    info = annotations.get(name)
    if not info: continue
    sp, sit, shared, iso = info
    col = sp_col(sp)
    y   = i + 0.5

    # alternating row tint
    if i % 2 == 0:
        ax_lab.axhspan(y-0.5, y+0.5, color="#F8F9FA", zorder=0)

    # bold + coloured isolate name
    ax_lab.text(0.01, y, name, color=col, fontsize=fsz, va="center",
                fontweight="bold" if (sit or shared) else "normal")

    # species in italic
    sp_display = sp if sp != "unknown" else "—"
    ax_lab.text(0.42, y, sp_display, color=col, fontsize=fsz,
                va="center", style="italic")

    # flags
    if sit and shared:
        ax_lab.text(0.80, y, "🔴 Sit  ★", color="#DC2626",
                    fontsize=fsz, va="center", fontweight="bold")
    elif sit:
        ax_lab.text(0.80, y, "🔴 Sit", color="#DC2626",
                    fontsize=fsz, va="center", fontweight="bold")
    elif shared:
        ax_lab.text(0.80, y, "★ shared", color="#DB2777",
                    fontsize=fsz, va="center", fontweight="bold")

ax_lab.set_title("Annotations", color=DARK, fontsize=11,
                 fontweight="bold", loc="left", pad=12)

# ── legend ────────────────────────────────────────────────────────────────────
present_sp = sorted({annotations[nm][0]
                     for nm in term_order if nm in annotations
                     and annotations[nm][0] != "unknown"})
leg = [mpatches.Patch(color=sp_col(s), label=s) for s in present_sp]
leg += [mpatches.Patch(color="#CBD5E1", label="unknown species"),
        mpatches.Patch(color="#DC2626", label="🔴  Sit iron transport genes"),
        mpatches.Patch(color="#DB2777", label="★   Shared between women")]
ax_tree.legend(handles=leg, fontsize=8, loc="lower right",
               frameon=True, edgecolor=BORDER, facecolor=BG, ncol=2,
               title="Host species & markers", title_fontsize=9,
               borderpad=1)

plt.tight_layout(pad=2)
plt.savefig(OUT_TREE, dpi=160, bbox_inches="tight", facecolor=BG)
print(f"  Saved → {OUT_TREE}")

# =============================================================================
# STEP 6 — DISTANCE TABLE
# =============================================================================
print("\nStep 6: Saving distances...")

def get_dist(a, b):
    try:    return dist_matrix[a, b]
    except Exception:
        try: return dist_matrix[b, a]
        except: return float("nan")

with open(OUT_DIST, "w") as fh:
    fh.write("isolate_a\tisolate_b\tspecies_a\tspecies_b\t"
             "distance\tboth_sit\tshared_phage\tinterpretation\n")
    for i in range(len(dm_ids)):
        for j in range(i+1, len(dm_ids)):
            a, b = dm_ids[i], dm_ids[j]
            ia = annotations.get(a, ("unknown",False,False,a))
            ib = annotations.get(b, ("unknown",False,False,b))
            d  = get_dist(a, b)
            interp = ("identical"  if d<0.01 else
                      "very close" if d<0.05 else
                      "related"    if d<0.15 else "distant")
            fh.write(f"{a}\t{b}\t{ia[0]}\t{ib[0]}\t{d:.5f}\t"
                     f"{'yes' if ia[1] and ib[1] else 'no'}\t"
                     f"{'yes' if ia[2] and ib[2] else 'no'}\t"
                     f"{interp}\n")

print(f"  Saved → {OUT_DIST}")

print("\n=== KEY DISTANCES ===")
print(f"  {'Pair':<50} {'Dist':>7}  Verdict")
print("  " + "-"*70)
for i in range(len(dm_ids)):
    for j in range(i+1, len(dm_ids)):
        a, b = dm_ids[i], dm_ids[j]
        ia = annotations.get(a,("?",False,False,a))
        ib = annotations.get(b,("?",False,False,b))
        d  = get_dist(a, b)
        label = ""
        if ia[1] and ib[1]:   label = "[Sit×Sit]"
        elif ia[2] and ib[2]: label = "[shared] "
        if label:
            interp = ("identical"  if d<0.01 else
                      "very close" if d<0.05 else
                      "related"    if d<0.15 else "distant")
            print(f"  {label}  {a:<22} vs {b:<22} {d:>7.4f}  {interp}")

print("\n=== DONE ===")
print(f"Tree  : {OUT_TREE}")
print(f"Table : {OUT_DIST}")
print()
print("Copy to Windows — run in WSL (not on Edzna):")
print("  scp edzna:~/utir/results/utir_terl_phylogeny_v3.png "
      "/mnt/c/Users/User/Desktop/")
print("  scp edzna:~/utir/results/utir_terl_distances_v3.tsv "
      "/mnt/c/Users/User/Desktop/")
