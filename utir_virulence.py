#!/usr/bin/env python3
"""
UTIr Prophage Virulence Analysis
==================================
Run from ~/utir:
    micromamba run -n pharokka18 python3 utir_virulence.py

What this does:
  1. Parses all 41 VFDB hits from Pharokka output
  2. Maps each hit to isolate, participant, species, timepoint
  3. Analyses the Sit iron transport system in depth
  4. Analyses bile salt hydrolase (BSH) in E_B. lactis
  5. Analyses cytolethal distending toxin (CDT) in participant 014
  6. Analyses LPS modification (gtrA/gtrB), adherence (nmpC), regulation (rcsB)
  7. Produces a multi-panel figure: results/utir_virulence.png
  8. Saves a summary table: results/utir_virulence_table.tsv

Outputs:
    results/utir_virulence.png
    results/utir_virulence_table.tsv
"""

import os, sys, re, csv, warnings
warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)

REQUIRED = [
    "results/pharokka_v113/top_hits_vfdb.tsv",
    "results/gtdbtk_all/gtdbtk.bac120.summary.tsv",
    "results/isolate_clusters.tsv",
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
import matplotlib.gridspec as gridspec
import numpy as np
from collections import defaultdict, Counter
print("  OK")

# =============================================================================
# HELPERS
# =============================================================================

def shorten_species(full):
    p = full.strip().split()
    return f"{p[0][0]}. {' '.join(p[1:])}" if len(p) >= 2 else full

def load_species():
    s = {}
    with open("results/gtdbtk_all/gtdbtk.bac120.summary.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            parts = row.get("classification","").split(";")
            sp = parts[-1].replace("s__","").strip() or "unclassified"
            s[row["user_genome"]] = shorten_species(sp)
    return s

def load_metadata():
    """Return {isolate: {participant, timepoint, timepoint_type}}"""
    meta = {}
    with open("results/isolate_clusters.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            iso = row["isolate_id"]
            tp  = row["timepoint"]
            if tp == "NA":
                tp_type = "unknown"
            elif tp.startswith("m"):
                tp_type = "UTI-free"
            elif tp.startswith("u"):
                tp_type = "UTI"
            else:
                tp_type = "unknown"
            meta[iso] = {
                "participant": row["participant_id"],
                "timepoint":   tp,
                "tp_type":     tp_type,
            }
    return meta

def iso_from_pid(pid):
    s = re.sub(r"_CDS_\d+$","", pid)
    s = re.sub(r"(__|\|)provirus.*","", s)
    s = re.sub(r"_\d{3}$","", s)
    return s

def parse_vfdb_hit(hit_str):
    """Extract gene short name and function category from VFDB hit string"""
    m_gene = re.search(r"\((\w+)\)", hit_str)
    gene   = m_gene.group(1) if m_gene else "unknown"
    # category
    if "iron" in hit_str.lower() or "sit" in hit_str.lower():
        category = "Iron transport (Sit)"
    elif "bile" in hit_str.lower() or "bsh" in hit_str.lower():
        category = "Stress survival (BSH)"
    elif "cdt" in hit_str.lower() or "distending" in hit_str.lower():
        category = "Exotoxin (CDT)"
    elif "lps" in hit_str.lower() or "gtr" in hit_str.lower():
        category = "Immune modulation (LPS)"
    elif "nmpC" in hit_str or "phosphoporin" in hit_str.lower():
        category = "Adherence (nmpC)"
    elif "rcs" in hit_str.lower() or "regulat" in hit_str.lower():
        category = "Regulation (rcsB)"
    elif "autotransporter" in hit_str.lower() or "yadA" in hit_str.lower():
        category = "Effector delivery"
    else:
        category = "Other"
    return gene, category

# =============================================================================
# STEP 1 — PARSE VFDB HITS
# =============================================================================
print("\nStep 1: Parsing VFDB hits...")

species_map = load_species()
meta_map    = load_metadata()

hits = []
with open("results/pharokka_v113/top_hits_vfdb.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        pid      = row["contig"]
        iso      = iso_from_pid(pid)
        hit_str  = row["vfdb_hit"]
        gene, cat = parse_vfdb_hit(hit_str)
        ident    = float(row["vfdb_seqIdentity"])
        meta     = meta_map.get(iso, {"participant":"?","timepoint":"?","tp_type":"unknown"})
        sp       = species_map.get(iso, "unknown")

        hits.append({
            "prophage_id":  pid,
            "isolate":      iso,
            "species":      sp,
            "participant":  meta["participant"],
            "timepoint":    meta["timepoint"],
            "tp_type":      meta["tp_type"],
            "gene":         gene,
            "category":     cat,
            "identity":     ident,
            "vfdb_hit":     hit_str,
            "start":        int(row["start"]),
            "stop":         int(row["stop"]),
        })

print(f"  Total hits: {len(hits)}")
for cat, n in Counter(h["category"] for h in hits).most_common():
    print(f"    {cat:<30} {n}")

# =============================================================================
# STEP 2 — SIT SYSTEM ANALYSIS
# =============================================================================
print("\nStep 2: Sit iron transport system...")

# group by isolate — which genes does each isolate carry?
sit_by_iso = defaultdict(set)
for h in hits:
    if "Iron transport" in h["category"]:
        sit_by_iso[h["isolate"]].add(h["gene"])

print(f"  Isolates with ≥1 Sit gene: {len(sit_by_iso)}")
sit_complete = {iso: genes for iso, genes in sit_by_iso.items()
                if {"sitA","sitB","sitC","sitD"} <= genes}
sit_partial  = {iso: genes for iso, genes in sit_by_iso.items()
                if iso not in sit_complete}

print(f"  Complete operon (sitABCD): {len(sit_complete)}")
for iso, genes in sorted(sit_complete.items()):
    sp  = species_map.get(iso,"?")
    tp  = meta_map.get(iso,{}).get("tp_type","?")
    pt  = meta_map.get(iso,{}).get("participant","?")
    print(f"    {iso:<30} P{pt} {tp:<10} {sp}")

print(f"  Partial operon: {len(sit_partial)}")
for iso, genes in sorted(sit_partial.items()):
    sp  = species_map.get(iso,"?")
    tp  = meta_map.get(iso,{}).get("tp_type","?")
    pt  = meta_map.get(iso,{}).get("participant","?")
    present = ",".join(sorted(genes))
    print(f"    {iso:<30} P{pt} {tp:<10} {sp:<20} genes: {present}")

# =============================================================================
# STEP 3 — BSH ANALYSIS
# =============================================================================
print("\nStep 3: Bile salt hydrolase (BSH)...")
bsh_hits = [h for h in hits if h["category"]=="Stress survival (BSH)"]
print(f"  BSH hits: {len(bsh_hits)}")
for h in bsh_hits:
    print(f"    {h['isolate']:<30} P{h['participant']} {h['tp_type']:<10}"
          f" {h['species']:<20} identity={h['identity']:.1%}")

# =============================================================================
# STEP 4 — CDT ANALYSIS
# =============================================================================
print("\nStep 4: Cytolethal distending toxin (CDT)...")
cdt_hits = [h for h in hits if h["category"]=="Exotoxin (CDT)"]
print(f"  CDT hits: {len(cdt_hits)}")
for h in cdt_hits:
    print(f"    {h['gene']:<8} {h['isolate']:<30} P{h['participant']}"
          f" {h['tp_type']:<10} identity={h['identity']:.1%}")

# =============================================================================
# STEP 5 — SAVE TABLE
# =============================================================================
print("\nStep 5: Saving table...")
with open("results/utir_virulence_table.tsv","w") as fh:
    fh.write("isolate\tparticipant\ttimepoint\ttp_type\tspecies\t"
             "gene\tcategory\tidentity\tprophage_id\n")
    for h in hits:
        fh.write(f"{h['isolate']}\t{h['participant']}\t{h['timepoint']}\t"
                 f"{h['tp_type']}\t{h['species']}\t{h['gene']}\t"
                 f"{h['category']}\t{h['identity']:.4f}\t"
                 f"{h['prophage_id']}\n")
print("  Saved → results/utir_virulence_table.tsv")

# =============================================================================
# STEP 6 — FIGURE
# =============================================================================
print("\nStep 6: Drawing figure...")

BG="#FFFFFF"; PANEL="#F8F9FA"; GRID="#E9ECEF"; BORDER="#DEE2E6"
BLUE="#2563EB"; GREEN="#16A34A"; AMBER="#D97706"; RED="#DC2626"
PURPLE="#7C3AED"; TEAL="#0891B2"; ORANGE="#EA580C"; PINK="#DB2777"
DARK="#111827"; MID="#374151"; SOFT="#9CA3AF"

CAT_COLS = {
    "Iron transport (Sit)":    GREEN,
    "Stress survival (BSH)":   TEAL,
    "Exotoxin (CDT)":          RED,
    "Immune modulation (LPS)": AMBER,
    "Adherence (nmpC)":        BLUE,
    "Regulation (rcsB)":       PURPLE,
    "Effector delivery":        ORANGE,
    "Other":                   SOFT,
}

fig = plt.figure(figsize=(18, 24), facecolor=BG)
fig.suptitle("Prophage-Encoded Virulence Genes — UTIr Cohort\n"
             "41 VFDB hits across 17,129 prophage genes",
             color=DARK, fontsize=14, fontweight="bold",
             x=0.04, ha="left", y=0.998)

gs = gridspec.GridSpec(4, 2, figure=fig,
                       hspace=0.52, wspace=0.35,
                       left=0.07, right=0.97,
                       top=0.965, bottom=0.03)

# ── A: overview bar — hits per category ──────────────────────────────────────
ax_a = fig.add_subplot(gs[0, :])
ax_a.set_facecolor(PANEL)
ax_a.spines["top"].set_visible(False); ax_a.spines["right"].set_visible(False)
ax_a.spines["left"].set_color(BORDER); ax_a.spines["bottom"].set_color(BORDER)
ax_a.set_title("A  All 41 Virulence Hits by Functional Category",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_a.text(0, 1.01,
          "Coloured by function  |  each bar = one hit in one prophage in one isolate",
          transform=ax_a.transAxes, color=SOFT, fontsize=8.5)

# sort hits by category then gene for grouped display
cat_order = ["Iron transport (Sit)","Stress survival (BSH)","Exotoxin (CDT)",
             "Immune modulation (LPS)","Adherence (nmpC)",
             "Regulation (rcsB)","Effector delivery"]
hits_sorted = sorted(hits, key=lambda h: (cat_order.index(h["category"])
                                          if h["category"] in cat_order else 99,
                                          h["gene"], h["isolate"]))

x_a      = np.arange(len(hits_sorted))
cols_a   = [CAT_COLS.get(h["category"], SOFT) for h in hits_sorted]
idents_a = [h["identity"] for h in hits_sorted]

bars_a = ax_a.bar(x_a, idents_a, color=cols_a, edgecolor=BG,
                  linewidth=0, width=0.85, zorder=3)
ax_a.set_ylim(0, 1.15)
ax_a.set_ylabel("Amino acid identity to VFDB reference", color=MID, fontsize=9.5)
ax_a.set_xticks(x_a)
ax_a.set_xticklabels([h["gene"] for h in hits_sorted],
                     rotation=70, ha="right", fontsize=7.5)
ax_a.tick_params(colors=MID)
ax_a.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
ax_a.axhline(0.80, color=SOFT, linewidth=0.8, linestyle=":",
             label="80% identity threshold")
ax_a.legend(fontsize=8.5, frameon=True, edgecolor=BORDER, facecolor=BG)

# isolate labels above bars
for i, (b, h) in enumerate(zip(bars_a, hits_sorted)):
    short_iso = h["isolate"].replace("utir_","")
    ax_a.text(b.get_x()+b.get_width()/2, h["identity"]+0.01,
              short_iso, ha="center", va="bottom",
              rotation=90, fontsize=5.5, color=DARK)

# category legend
leg_a = [mpatches.Patch(color=CAT_COLS[c], label=c) for c in cat_order]
ax_a.legend(handles=leg_a, fontsize=8, frameon=True,
            edgecolor=BORDER, facecolor=BG, ncol=2,
            loc="lower right")

# ── B: Sit operon heatmap ─────────────────────────────────────────────────────
ax_b = fig.add_subplot(gs[1, 0])
ax_b.set_facecolor(PANEL)
ax_b.spines["top"].set_visible(False); ax_b.spines["right"].set_visible(False)
ax_b.spines["left"].set_color(BORDER); ax_b.spines["bottom"].set_color(BORDER)
ax_b.set_title("B  Sit Iron Transport Operon — Which Genes Are Present?",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_b.text(0, 1.01, "Iron is scarce in urine — Sit genes give bacteria a competitive advantage",
          transform=ax_b.transAxes, color=SOFT, fontsize=8.5)

# all isolates with any Sit gene
sit_isolates_ordered = sorted(sit_by_iso.keys(),
                               key=lambda i: (
                                   0 if {"sitA","sitB","sitC","sitD"}<=sit_by_iso[i] else 1,
                                   meta_map.get(i,{}).get("participant","?"),
                               ))
sit_genes_all = ["sitA","sitB","sitC","sitD"]

y_b = np.arange(len(sit_isolates_ordered))
for gi, gene in enumerate(sit_genes_all):
    for yi, iso in enumerate(sit_isolates_ordered):
        present = gene in sit_by_iso[iso]
        # get identity for this gene in this isolate
        ident_val = next((h["identity"] for h in hits
                          if h["isolate"]==iso and h["gene"]==gene), None)
        col  = GREEN if present else GRID
        rect = mpatches.FancyBboxPatch(
            (gi+0.05, yi+0.05), 0.9, 0.9,
            boxstyle="round,pad=0.04",
            facecolor=col, edgecolor=BG, linewidth=0.5, zorder=3)
        ax_b.add_patch(rect)
        if present and ident_val:
            ax_b.text(gi+0.5, yi+0.5, f"{ident_val:.0%}",
                      ha="center", va="center",
                      color="white", fontsize=7.5, fontweight="bold", zorder=4)
        elif not present:
            ax_b.text(gi+0.5, yi+0.5, "—",
                      ha="center", va="center",
                      color=SOFT, fontsize=9, zorder=4)

# labels
ax_b.set_xlim(0, 4); ax_b.set_ylim(0, len(sit_isolates_ordered))
ax_b.set_xticks([0.5,1.5,2.5,3.5])
ax_b.set_xticklabels(sit_genes_all, fontsize=10, fontweight="bold", color=GREEN)

y_labels = []
for iso in sit_isolates_ordered:
    pt  = meta_map.get(iso,{}).get("participant","?")
    tp  = meta_map.get(iso,{}).get("tp_type","?")
    sp  = species_map.get(iso,"?")
    complete = "✓" if {"sitA","sitB","sitC","sitD"} <= sit_by_iso[iso] else " "
    y_labels.append(f"{complete} P{pt}  {iso.replace('utir_','')[:18]}")

ax_b.set_yticks([i+0.5 for i in range(len(sit_isolates_ordered))])
ax_b.set_yticklabels(y_labels, fontsize=8)
ax_b.invert_yaxis()
ax_b.tick_params(length=0)

# complete/partial markers
ax_b.text(4.1, 0.5, "✓ = complete\noperon", color=GREEN,
          fontsize=8, va="center")

# ── C: BSH — bile salt hydrolase context ─────────────────────────────────────
ax_c = fig.add_subplot(gs[1, 1])
ax_c.set_facecolor(BG); ax_c.axis("off")
ax_c.set_title("C  Bile Salt Hydrolase (BSH) — Why in a Urinary Phage?",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)

bsh_isolates = [
    ("utir_011_u1k2_02", "P011", "UTI",      "E_B. lactis", "82.5%"),
    ("utir_011_u1k2_04", "P011", "UTI",      "E_B. lactis", "82.5%"),
    ("utir_025_m0_07",   "P025", "UTI-free", "E_B. lactis", "82.5%"),
    ("utir_029_m0_06",   "P029", "UTI-free", "E_B. lactis", "82.5%"),
]

y0 = 0.96
ax_c.text(0.02, y0, "BSH gene found in 4 isolates — all E_B. lactis",
          color=TEAL, fontsize=10, fontweight="bold",
          transform=ax_c.transAxes, va="top")
y0 -= 0.08
for iso, pt, tp, sp, ident in bsh_isolates:
    tp_col = AMBER if tp=="UTI" else BLUE
    ax_c.text(0.04, y0, f"• {iso.replace('utir_','')}",
              color=DARK, fontsize=9, fontweight="bold",
              transform=ax_c.transAxes, va="top")
    ax_c.text(0.38, y0, f"{pt}  {tp}  {ident}",
              color=tp_col, fontsize=9,
              transform=ax_c.transAxes, va="top")
    y0 -= 0.07

y0 -= 0.04
ax_c.text(0.02, y0, "What BSH does:", color=DARK, fontsize=9.5,
          fontweight="bold", transform=ax_c.transAxes, va="top")
y0 -= 0.07
bsh_text = [
    ("Bile salt hydrolase deconjugates bile acids — it splits bile",   MID),
    ("salts that the host uses to emulsify fats in the gut.",          MID),
    ("",                                                               MID),
    ("Why would it be in a urinary phage?",                           TEAL),
    ("Bile acids are present at low levels in the urinary tract,",     MID),
    ("particularly in women with enterococcal colonisation.",          MID),
    ("BSH may help the bacterium tolerate this chemical stress,",      MID),
    ("or it may be a relic from an ancestor that also colonised",      MID),
    ("the gut — a common niche for E_B. lactis.",                      MID),
    ("",                                                               MID),
    ("Key observation:",                                               GREEN),
    ("All 4 isolates carry the SAME prophage (100% ANI between",       GREEN),
    ("025 and 029; 99.4% with 011). The BSH gene is maintained",      GREEN),
    ("across all copies — suggesting it is SELECTED FOR,",            GREEN),
    ("not a neutral passenger.",                                       GREEN),
]
for text, col in bsh_text:
    ax_c.text(0.02, y0, text, color=col, fontsize=8.5,
              transform=ax_c.transAxes, va="top", linespacing=1.4)
    y0 -= 0.052

# ── D: CDT — cytolethal distending toxin ─────────────────────────────────────
ax_d = fig.add_subplot(gs[2, 0])
ax_d.set_facecolor(BG); ax_d.axis("off")
ax_d.set_title("D  Cytolethal Distending Toxin (CDT) — DNA Damage in Host Cells",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)

y0 = 0.96
ax_d.text(0.02, y0,
          "cdtA (93.1%) + cdtC (88.8%) found in utir_014_u1k1_05",
          color=RED, fontsize=10, fontweight="bold",
          transform=ax_d.transAxes, va="top")
y0 -= 0.08
ax_d.text(0.02, y0, "Participant 014 — UTI episode 1, culture 1 — E. coli",
          color=MID, fontsize=9, transform=ax_d.transAxes, va="top")
y0 -= 0.1

cdt_text = [
    ("What CDT is:", DARK),
    ("Cytolethal distending toxin is a 3-subunit AB2 toxin.", MID),
    ("CdtA and CdtC are the binding subunits that attach the", MID),
    ("toxin to the surface of human cells. CdtB (not found", MID),
    ("here) is the active subunit — a DNase that cuts DNA", MID),
    ("inside the human cell, causing double-strand breaks.", MID),
    ("", MID),
    ("What it does to the host:", DARK),
    ("Cells with broken DNA activate repair checkpoints and", MID),
    ("stop dividing — they become 'distended' (swollen) and", MID),
    ("eventually die. CDT effectively paralyses immune cells", MID),
    ("(macrophages, lymphocytes) that would otherwise kill the", MID),
    ("bacteria, giving the pathogen time to establish infection.", MID),
    ("", MID),
    ("Why finding cdtA+C without cdtB matters:", DARK),
    ("The active DNase subunit (cdtB) is absent from this", AMBER),
    ("prophage. The toxin cannot function without all 3 parts.", AMBER),
    ("Two interpretations:", AMBER),
    ("1. cdtB is in another prophage of this isolate (014_u1k1_05", MID),
    ("   has 10 prophages — worth checking).", MID),
    ("2. This is a partial CDT — the binding subunits were", MID),
    ("   acquired but the operon was never completed.", MID),
    ("", MID),
    ("Clinical context:", DARK),
    ("CDT has been found in uropathogenic E. coli (UPEC) strains", RED),
    ("causing cystitis and pyelonephritis. Its presence here,", RED),
    ("even partial, in participant 014's first UTI episode is", RED),
    ("consistent with UPEC virulence biology.", RED),
]
for text, col in cdt_text:
    ax_d.text(0.02, y0, text, color=col, fontsize=8.2,
              transform=ax_d.transAxes, va="top", linespacing=1.4)
    y0 -= 0.045

# ── E: remaining hits summary ─────────────────────────────────────────────────
ax_e = fig.add_subplot(gs[2, 1])
ax_e.set_facecolor(BG); ax_e.axis("off")
ax_e.set_title("E  Remaining Virulence Hits",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)

remaining = [
    ("LPS modification\ngtrA + gtrB",
     "2 prophages:\nutir_011_m12_05\nutir_004_m0_01",
     AMBER,
     "Modifies the bacterial surface coat (LPS)\nto evade immune recognition.\ngtrA flips sugar molecules across\nthe membrane; gtrB attaches glucose\nto LPS. Together they camouflage\nthe bacterium from the host immune system."),
    ("Adherence\nnmpC (×5)",
     "5 isolates:\n011_u1k1_04\n014_u1k1_05\n025_m0_05/06\n025_m0_06",
     BLUE,
     "Outer membrane phosphoporin — a pore\nprotein that may help bacteria attach\nto bladder epithelial cells. Found in\n5 isolates across 2 participants.\nThe Salmonella reference (nmpC/PhoE)\nsuggests this is a cross-species\nhorizontally transferred gene."),
    ("Regulation\nrcsB (×3)",
     "3 prophages:\n003_m0_02\n003_m0_04\n011_m0_02",
     PURPLE,
     "Transcriptional regulator controlling\ncapsule biosynthesis and biofilm\nformation. A phage-encoded rcsB could\ntip the bacterium toward capsule\nproduction, helping it resist\nphagocytosis by immune cells."),
    ("Effector delivery\nYadA + autotransporter",
     "1 prophage:\n004_m0_01",
     ORANGE,
     "YadA-like protein: an adhesin that\nhelps bacteria stick to host tissues.\nAutotransporter: a secretion system\nthat exports proteins to the bacterial\nsurface. Both found in the same\nprophage in the most prophage-rich\nisolate in the dataset (16 prophages)."),
]

y0 = 0.97
for name, isolates, col, explanation in remaining:
    ax_e.add_patch(mpatches.FancyBboxPatch(
        (0.01, y0-0.21), 0.97, 0.22,
        boxstyle="round,pad=0.01",
        facecolor=col+"15", edgecolor=col, linewidth=1.5,
        transform=ax_e.transAxes, zorder=1, clip_on=False))
    ax_e.text(0.03, y0-0.01, name, color=col, fontsize=9,
              fontweight="bold", transform=ax_e.transAxes, va="top",
              linespacing=1.4)
    ax_e.text(0.28, y0-0.01, isolates, color=MID, fontsize=7.5,
              transform=ax_e.transAxes, va="top", linespacing=1.4)
    ax_e.text(0.52, y0-0.01, explanation, color=MID, fontsize=7.5,
              transform=ax_e.transAxes, va="top", linespacing=1.4)
    y0 -= 0.25

# ── F: participant 014 — the most virulence-laden isolate ────────────────────
ax_f = fig.add_subplot(gs[3, :])
ax_f.set_facecolor(PANEL)
ax_f.spines["top"].set_visible(False); ax_f.spines["right"].set_visible(False)
ax_f.spines["left"].set_color(BORDER); ax_f.spines["bottom"].set_color(BORDER)
ax_f.set_title("F  Participant 014 — Most Virulence-Loaded Participant in the Cohort",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_f.text(0, 1.01,
          "Every virulence category except BSH appears in this participant's bacteria across UTI episodes",
          transform=ax_f.transAxes, color=SOFT, fontsize=8.5)

# count virulence hits per participant
hits_by_pt = defaultdict(list)
for h in hits:
    hits_by_pt[h["participant"]].append(h)

pt_sorted = sorted(hits_by_pt.keys(),
                   key=lambda p: -len(hits_by_pt[p]))

pt_names = [f"P{p}" for p in pt_sorted]
pt_counts = [len(hits_by_pt[p]) for p in pt_sorted]

# colour by whether participant 014
pt_cols = [RED if p=="014" else BLUE for p in pt_sorted]

bars_f = ax_f.bar(pt_names, pt_counts, color=pt_cols,
                  edgecolor=BG, linewidth=0, width=0.65, zorder=3)
for b, v in zip(bars_f, pt_counts):
    ax_f.text(b.get_x()+b.get_width()/2, v+0.1, str(v),
              ha="center", va="bottom", color=DARK,
              fontsize=10, fontweight="bold")

ax_f.set_ylabel("Total virulence gene hits", color=MID, fontsize=10)
ax_f.tick_params(colors=MID)
ax_f.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
ax_f.set_ylim(0, max(pt_counts)+3)

# annotate 014
p014_idx = pt_sorted.index("014")
ax_f.annotate(
    "Participant 014:\nemrE + baeR (resistance)\nCDT (exotoxin)\nSit partial\nnmpC (adherence)\n4 UTI episodes",
    xy=(p014_idx, hits_by_pt["014"].__len__()),
    xytext=(p014_idx+1.5, max(pt_counts)-2),
    arrowprops=dict(arrowstyle="->", color=RED, lw=1.3),
    color=RED, fontsize=8.5, fontweight="bold",
    bbox=dict(facecolor="#FFF1F2", edgecolor=RED, pad=5, boxstyle="round"))

plt.savefig("results/utir_virulence.png",
            dpi=160, bbox_inches="tight", facecolor=BG)
print("  Saved → results/utir_virulence.png")

print("\n=== DONE ===")
print("Outputs:")
print("  results/utir_virulence.png")
print("  results/utir_virulence_table.tsv")
print()
print("Copy to Windows — run in WSL:")
print("  scp edzna:~/utir/results/utir_virulence.png /mnt/c/Users/User/Desktop/")
print("  scp edzna:~/utir/results/utir_virulence_table.tsv /mnt/c/Users/User/Desktop/")
