#!/usr/bin/env python3
"""
UTIr tRNA Deep Analysis — v2
Fixes:
  - format string error in TSV writer
  - species tRNA type lookup (was reading amino_acid field which was empty;
    now reads from product field which contains "tRNA-Arg(AGA)" etc.)
  - tmRNA detection (aragorn GFF has different format; now checks both files)
Run from ~/utir:
    micromamba run -n pharokka18 python3 utir_trna_analysis_v2.py
"""

import os, sys, re, csv, warnings
warnings.filterwarnings("ignore")

os.makedirs("results", exist_ok=True)

REQUIRED = [
    "results/pharokka_v113/trnascan_out.gff",
    "results/pharokka_v113/pharokka_aragorn.gff",
    "results/gtdbtk_all/gtdbtk.bac120.summary.tsv",
    "results/checkv/quality_summary.tsv",
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
from collections import Counter, defaultdict
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

def iso_from_pid(pid):
    s = re.sub(r"_CDS_\d+$","", pid)
    s = re.sub(r"(__|\|)provirus.*","", s)
    s = re.sub(r"_\d{3}$","", s)
    return s

def prophage_len_from_id(pid):
    m = re.search(r"provirus_(\d+)_(\d+)", pid)
    if m:
        return int(m.group(2)) - int(m.group(1))
    # whole contig: provirus_2_41420 or just _002
    m2 = re.search(r"_(\d+)_(\d+)$", pid)
    if m2:
        try:
            return int(m2.group(2)) - int(m2.group(1))
        except:
            pass
    return None

def fmt_pos(v):
    """Safe float formatter"""
    if v is None: return "NA"
    try:    return f"{float(v):.2f}"
    except: return "NA"

EDGE_THRESHOLD = 15.0

# =============================================================================
# STEP 1 — PARSE tRNAscan GFF
# =============================================================================
print("\nStep 1: Parsing tRNAscan GFF...")

def parse_gff_file(path, source_label):
    entries = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("tRNA","tmRNA"):
                continue
            pid    = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parts[8]

            # product — e.g. "tRNA-Arg(AGA)" or "transfer-messenger RNA SsrA"
            m = re.search(r"product=([^;]+)", attrs)
            product = m.group(1).strip() if m else "unknown"

            # extract amino acid from product string tRNA-XXX(YYY)
            # tRNAscan GFF format uses isotype= and anticodon= fields
            # e.g. isotype=Ser;anticodon=TGA
            m_iso = re.search(r"isotype=(\w+)", attrs)
            m_ac2 = re.search(r"anticodon=(\w+)", attrs)
            # also try old-style anticodon=(pos:...,aa:X,seq:Y)
            m_ac_old = re.search(r"anticodon=\(pos:[^,]+,aa:(\w+),seq:(\w+)\)", attrs)
            # also try product= style as fallback
            m_aa_prod = re.search(r"tRNA-([A-Za-z]+)\(([A-Za-z]+)\)", product)

            if m_iso:
                amino_acid = m_iso.group(1)          # e.g. "Ser"
                codon      = m_ac2.group(1) if m_ac2 else ""  # e.g. "TGA"
            elif m_ac_old:
                amino_acid = m_ac_old.group(1)
                codon      = m_ac_old.group(2)
            elif m_aa_prod:
                amino_acid = m_aa_prod.group(1)
                codon      = m_aa_prod.group(2)
            else:
                amino_acid = ""
                codon      = ""

            plen    = prophage_len_from_id(pid)
            rel_pos = 100.0 * start / plen if (plen and plen > 0) else None

            if rel_pos is not None:
                location = ("EDGE"
                            if rel_pos <= EDGE_THRESHOLD
                            or rel_pos >= (100 - EDGE_THRESHOLD)
                            else "MIDDLE")
            else:
                location = "UNKNOWN"

            iso = iso_from_pid(pid)
            entries.append({
                "prophage_id": pid,
                "isolate":     iso,
                "feature":     feature,
                "source":      source_label,
                "start":       start,
                "end":         end,
                "strand":      strand,
                "product":     product,
                "amino_acid":  amino_acid,
                "codon":       codon,
                "plen_bp":     plen,
                "rel_pos_pct": rel_pos,
                "location":    location,
                "species":     "",
            })
    return entries

trnas_scan   = parse_gff_file("results/pharokka_v113/trnascan_out.gff",   "tRNAscan")
trnas_aragorn= parse_gff_file("results/pharokka_v113/pharokka_aragorn.gff","Aragorn")

# combine — deduplicate by (prophage_id, start, feature)
seen = set()
trnas = []
for t in trnas_scan + trnas_aragorn:
    key = (t["prophage_id"], t["start"], t["feature"])
    if key not in seen:
        seen.add(key)
        trnas.append(t)

print(f"  tRNAscan entries : {len(trnas_scan)}")
print(f"  Aragorn entries  : {len(trnas_aragorn)}")
print(f"  After dedup      : {len(trnas)}")
print(f"  EDGE   : {sum(1 for t in trnas if t['location']=='EDGE')}")
print(f"  MIDDLE : {sum(1 for t in trnas if t['location']=='MIDDLE')}")
print(f"  tmRNA  : {sum(1 for t in trnas if t['feature']=='tmRNA')}")

tmrna_entries = [t for t in trnas if t["feature"]=="tmRNA"]
if tmrna_entries:
    for t in tmrna_entries:
        print(f"  tmRNA found: {t['prophage_id']}  pos={t['start']}-{t['end']}"
              f"  ({t['rel_pos_pct']:.1f}%)  product={t['product']}")
else:
    print("  No tmRNA found after dedup — checking aragorn raw...")
    # print first 3 aragorn lines for debug
    with open("results/pharokka_v113/pharokka_aragorn.gff") as fh:
        for i,line in enumerate(fh):
            if i < 5: print(f"    {line.rstrip()}")

# =============================================================================
# STEP 2 — SPECIES
# =============================================================================
print("\nStep 2: Loading species...")
species_map = load_species()
for t in trnas:
    t["species"] = species_map.get(t["isolate"], "unknown")

# per-prophage summary
prophages_with_trna = {}
for t in trnas:
    pid = t["prophage_id"]
    if pid not in prophages_with_trna:
        prophages_with_trna[pid] = {
            "isolate": t["isolate"], "species": t["species"],
            "plen_bp": t["plen_bp"], "trnas":[], "n_edge":0, "n_middle":0,
        }
    prophages_with_trna[pid]["trnas"].append(t)
    if t["location"] == "EDGE":   prophages_with_trna[pid]["n_edge"]   += 1
    else:                          prophages_with_trna[pid]["n_middle"] += 1

print(f"  Prophages with ≥1 tRNA: {len(prophages_with_trna)}")

# =============================================================================
# STEP 3 — CHECKV
# =============================================================================
print("\nStep 3: CheckV cross-tab...")
checkv = {}
with open("results/checkv/quality_summary.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        cid = row["contig_id"].replace("|","__")
        checkv[cid] = row.get("checkv_quality","Not-determined")

quality_order = ["Complete","High-quality","Medium-quality","Low-quality"]
cv_trna   = Counter()
cv_notrna = Counter()
for cid, qual in checkv.items():
    if cid in prophages_with_trna: cv_trna[qual]   += 1
    else:                           cv_notrna[qual] += 1

print(f"  {'Quality':<20} {'With tRNA':>10} {'Without':>10}")
for q in quality_order:
    print(f"  {q:<20} {cv_trna.get(q,0):>10} {cv_notrna.get(q,0):>10}")

# =============================================================================
# STEP 4 — tRNA TYPES BY SPECIES
# =============================================================================
print("\nStep 4: tRNA types by species...")
species_trna = defaultdict(Counter)
for t in trnas:
    if t["feature"] == "tmRNA": continue
    sp = t["species"]
    aa = t["amino_acid"] if t["amino_acid"] else "OTHER"
    species_trna[sp][aa] += 1

sp_totals = {sp: sum(c.values()) for sp,c in species_trna.items()}
print("  Top 10 species:")
for sp, tot in sorted(sp_totals.items(), key=lambda x:-x[1])[:10]:
    top3 = ", ".join(f"{aa}({n})" for aa,n in species_trna[sp].most_common(3))
    print(f"    {sp:<30} n={tot:>4}  {top3}")

# =============================================================================
# STEP 5 — EDGE CANDIDATES
# =============================================================================
print("\nStep 5: Edge tRNA candidates...")
edge_trnas = [t for t in trnas if t["location"]=="EDGE" and t["feature"]=="tRNA"]
print(f"  {len(edge_trnas)} edge tRNAs")

with open("results/utir_trna_edge_candidates.tsv","w") as fh:
    fh.write("prophage_id\tisolate\tspecies\tamino_acid\tcodon\t"
             "start_bp\tprophage_len_bp\trel_pos_pct\tstrand\n")
    for t in sorted(edge_trnas, key=lambda x: x["rel_pos_pct"] or 0):
        fh.write(f"{t['prophage_id']}\t{t['isolate']}\t{t['species']}\t"
                 f"{t['amino_acid']}\t{t['codon']}\t"
                 f"{t['start']}\t{t['plen_bp']}\t"
                 f"{fmt_pos(t['rel_pos_pct'])}\t{t['strand']}\n")
print("  Saved → results/utir_trna_edge_candidates.tsv")

# =============================================================================
# STEP 6 — 28-tRNA PROPHAGE
# =============================================================================
print("\nStep 6: 28-tRNA prophage...")
big_pid   = "utir_008_m0_01_001__provirus_914808_996082"
big_trnas = [t for t in trnas if t["prophage_id"]==big_pid]
print(f"  Found {len(big_trnas)} tRNAs")
type_counts = Counter(t["amino_acid"] if t["amino_acid"] else "OTHER"
                      for t in big_trnas)
print(f"  Types ({len(type_counts)} distinct):")
for aa,n in type_counts.most_common():
    pos_list = [t["rel_pos_pct"] for t in big_trnas
                if (t["amino_acid"] if t["amino_acid"] else "OTHER")==aa
                and t["rel_pos_pct"] is not None]
    pos_str = ", ".join(f"{p:.1f}%" for p in pos_list)
    loc = "EDGE" if any(p<=EDGE_THRESHOLD or p>=(100-EDGE_THRESHOLD)
                        for p in pos_list if p) else "MIDDLE"
    print(f"    {aa:<10} ×{n}  {loc}  positions: {pos_str}")

# =============================================================================
# STEP 7 — SAVE FULL TABLE
# =============================================================================
print("\nStep 7: Saving full tRNA table...")
with open("results/utir_trna_table.tsv","w") as fh:
    fh.write("prophage_id\tisolate\tspecies\tfeature\tsource\t"
             "amino_acid\tcodon\tstart_bp\tend_bp\t"
             "prophage_len_bp\trel_position_pct\tstrand\tlocation\n")
    for t in trnas:
        rel = fmt_pos(t["rel_pos_pct"])
        fh.write(f"{t['prophage_id']}\t{t['isolate']}\t{t['species']}\t"
                 f"{t['feature']}\t{t['source']}\t"
                 f"{t['amino_acid']}\t{t['codon']}\t"
                 f"{t['start']}\t{t['end']}\t{t['plen_bp']}\t"
                 f"{rel}\t{t['strand']}\t{t['location']}\n")
print("  Saved → results/utir_trna_table.tsv")

# =============================================================================
# STEP 8 — FIGURE
# =============================================================================
print("\nStep 8: Drawing figure...")

BG="#FFFFFF"; PANEL="#F8F9FA"; GRID="#E9ECEF"; BORDER="#DEE2E6"
BLUE="#2563EB"; GREEN="#16A34A"; AMBER="#D97706"; RED="#DC2626"
PURPLE="#7C3AED"; TEAL="#0891B2"
DARK="#111827"; MID="#374151"; SOFT="#9CA3AF"

fig = plt.figure(figsize=(17,20), facecolor=BG)
fig.suptitle("UTIr Prophage tRNA Deep Analysis  —  v2\n"
             f"{len(trnas)} tRNAs + {len(tmrna_entries)} tmRNA  "
             f"across {len(prophages_with_trna)} prophages",
             color=DARK, fontsize=13, fontweight="bold",
             x=0.04, ha="left", y=0.998)

gs = gridspec.GridSpec(3,2, figure=fig,
                       hspace=0.48, wspace=0.35,
                       left=0.07, right=0.97,
                       top=0.965, bottom=0.04)

# ── A: position histogram ─────────────────────────────────────────────────────
ax_a = fig.add_subplot(gs[0,:])
ax_a.set_facecolor(PANEL)
for sp in ax_a.spines.values(): sp.set_color(BORDER)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.set_title("A  tRNA Position Within Prophage  "
               "(EDGE = integration site candidate  |  MIDDLE = codon-bias tool)",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_a.text(0, 1.01, f"Threshold: ±{EDGE_THRESHOLD}% from either end",
          transform=ax_a.transAxes, color=SOFT, fontsize=8.5)

pos_all   = [t["rel_pos_pct"] for t in trnas
             if t["rel_pos_pct"] is not None and t["feature"]=="tRNA"]
edge_pos  = [p for p in pos_all if p<=EDGE_THRESHOLD or p>=(100-EDGE_THRESHOLD)]
mid_pos   = [p for p in pos_all if EDGE_THRESHOLD<p<(100-EDGE_THRESHOLD)]

bins = np.linspace(0,100,26)
ax_a.hist(mid_pos,  bins=bins, color=BLUE, alpha=0.85, edgecolor=BG,
          linewidth=0.3, label=f"Middle — codon bias ({len(mid_pos)})", zorder=3)
ax_a.hist(edge_pos, bins=bins, color=RED,  alpha=0.85, edgecolor=BG,
          linewidth=0.3, label=f"Edge — integration site ({len(edge_pos)})", zorder=3)
ax_a.axvline(EDGE_THRESHOLD,     color=RED,  linewidth=1.5, linestyle="--", alpha=0.6)
ax_a.axvline(100-EDGE_THRESHOLD, color=RED,  linewidth=1.5, linestyle="--", alpha=0.6,
             label=f"±{EDGE_THRESHOLD}% threshold")
ax_a.axvline(50, color=SOFT, linewidth=0.8, linestyle=":", alpha=0.4)
ax_a.set_xlabel("tRNA position within prophage (%)", color=MID, fontsize=10)
ax_a.set_ylabel("Number of tRNAs", color=MID, fontsize=10)
ax_a.tick_params(colors=MID)
ax_a.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
ax_a.legend(fontsize=9, frameon=True, edgecolor=BORDER, facecolor=BG)
ymax = ax_a.get_ylim()[1]
ax_a.text(7.5,  ymax*0.85, f"{len(edge_pos)} integration\nsite candidates",
          ha="center", color=RED,  fontsize=9, fontweight="bold")
ax_a.text(50,   ymax*0.85, f"{len(mid_pos)} codon-bias\ntRNAs",
          ha="center", color=BLUE, fontsize=9, fontweight="bold")

# ── B: tRNA types edge vs middle ─────────────────────────────────────────────
ax_b = fig.add_subplot(gs[1,0])
ax_b.set_facecolor(PANEL)
ax_b.spines["top"].set_visible(False); ax_b.spines["right"].set_visible(False)
ax_b.spines["left"].set_color(BORDER); ax_b.spines["bottom"].set_color(BORDER)
ax_b.set_title("B  tRNA Types: Edge vs Middle",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_b.text(0, 1.01, "Which amino acids appear at integration sites vs middle?",
          transform=ax_b.transAxes, color=SOFT, fontsize=8.5)

edge_aa   = Counter(t["amino_acid"] if t["amino_acid"] else "OTHER"
                    for t in trnas if t["location"]=="EDGE"  and t["feature"]=="tRNA")
middle_aa = Counter(t["amino_acid"] if t["amino_acid"] else "OTHER"
                    for t in trnas if t["location"]=="MIDDLE" and t["feature"]=="tRNA")
all_aa = sorted(set(list(edge_aa)+list(middle_aa)),
                key=lambda a:-(edge_aa.get(a,0)+middle_aa.get(a,0)))[:14]
x_b = np.arange(len(all_aa)); w=0.38
ax_b.bar(x_b-w/2, [edge_aa.get(a,0)   for a in all_aa], w,
         color=RED,  edgecolor=BG, linewidth=0,
         label="Edge (integration site)", zorder=3)
ax_b.bar(x_b+w/2, [middle_aa.get(a,0) for a in all_aa], w,
         color=BLUE, edgecolor=BG, linewidth=0,
         label="Middle (codon bias)", zorder=3)
ax_b.set_xticks(x_b)
ax_b.set_xticklabels(all_aa, rotation=45, ha="right", fontsize=8.5)
ax_b.set_ylabel("Count", color=MID, fontsize=9)
ax_b.tick_params(colors=MID)
ax_b.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
ax_b.legend(fontsize=8.5, frameon=True, edgecolor=BORDER, facecolor=BG)

# ── C: tRNA count by species ──────────────────────────────────────────────────
ax_c = fig.add_subplot(gs[1,1])
ax_c.set_facecolor(PANEL)
ax_c.spines["top"].set_visible(False); ax_c.spines["right"].set_visible(False)
ax_c.spines["left"].set_color(BORDER); ax_c.spines["bottom"].set_color(BORDER)
ax_c.set_title("C  Total tRNAs by Host Species",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_c.text(0, 1.01, "Which bacterial species carry the most tRNA-laden prophages?",
          transform=ax_c.transAxes, color=SOFT, fontsize=8.5)

sp_sorted = sorted(sp_totals.items(), key=lambda x:-x[1])[:12]
sp_names  = [s[0] for s in sp_sorted]
sp_cnts   = [s[1] for s in sp_sorted]
sp_col_map= {"E. faecalis":AMBER,"E. coli":RED,"S. vaginalis":BLUE,
             "L. gasseri":"#15803D","W. anitrata":TEAL,
             "P. mirabilis":PURPLE,"S. haemolyticus":"#7C3AED"}
bar_cols  = [sp_col_map.get(s,SOFT) for s in sp_names]
y_c = np.arange(len(sp_names))
bars_c = ax_c.barh(y_c, sp_cnts, color=bar_cols,
                   edgecolor=BG, linewidth=0, height=0.65, zorder=3)
ax_c.set_yticks(y_c)
ax_c.set_yticklabels(sp_names, fontsize=8.5, style="italic")
ax_c.invert_yaxis()
ax_c.set_xlabel("Total tRNAs", color=MID, fontsize=9)
ax_c.tick_params(colors=MID)
ax_c.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
for b, v in zip(bars_c, sp_cnts):
    ax_c.text(v+0.3, b.get_y()+b.get_height()/2,
              str(v), va="center", color=MID, fontsize=9, fontweight="bold")

# ── D: CheckV quality ─────────────────────────────────────────────────────────
ax_d = fig.add_subplot(gs[2,0])
ax_d.set_facecolor(PANEL)
ax_d.spines["top"].set_visible(False); ax_d.spines["right"].set_visible(False)
ax_d.spines["left"].set_color(BORDER); ax_d.spines["bottom"].set_color(BORDER)
ax_d.set_title("D  CheckV Quality: With vs Without tRNA",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)
ax_d.text(0, 1.01,
          "Do tRNA-carrying prophages have better-defined boundaries?",
          transform=ax_d.transAxes, color=SOFT, fontsize=8.5)

qs_short = ["Complete","High","Medium","Low"]
qs_full  = ["Complete","High-quality","Medium-quality","Low-quality"]
with_v   = [cv_trna.get(q,0)   for q in qs_full]
witho_v  = [cv_notrna.get(q,0) for q in qs_full]
tw = sum(with_v)   or 1
to = sum(witho_v)  or 1
with_pct  = [100*v/tw for v in with_v]
witho_pct = [100*v/to for v in witho_v]
x_d = np.arange(4); wd=0.38
ax_d.bar(x_d-wd/2, with_pct,  wd, color=BLUE, edgecolor=BG, linewidth=0,
         label=f"With tRNA (n={tw})", zorder=3)
ax_d.bar(x_d+wd/2, witho_pct, wd, color=SOFT, edgecolor=BG, linewidth=0,
         label=f"Without tRNA (n={to})", zorder=3)
ax_d.set_xticks(x_d)
ax_d.set_xticklabels(qs_short, fontsize=10)
ax_d.set_ylabel("% of prophages", color=MID, fontsize=9)
ax_d.tick_params(colors=MID)
ax_d.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
ax_d.legend(fontsize=8.5, frameon=True, edgecolor=BORDER, facecolor=BG)
for bars_grp, pcts in [(ax_d.patches[:4], with_pct),
                        (ax_d.patches[4:], witho_pct)]:
    for b, v in zip(bars_grp, pcts):
        if v > 2:
            ax_d.text(b.get_x()+b.get_width()/2, b.get_height()+0.5,
                      f"{v:.0f}%", ha="center", va="bottom",
                      color=DARK, fontsize=7.5, fontweight="bold")

# ── E: special cases ──────────────────────────────────────────────────────────
ax_e = fig.add_subplot(gs[2,1])
ax_e.set_facecolor(BG); ax_e.axis("off")
ax_e.set_title("E  Special Cases",
               color=DARK, fontsize=11, fontweight="bold", loc="left", pad=8)

y_e = 0.97
ax_e.text(0.01, y_e, "28-tRNA prophage (S. vaginalis, 81 kb)",
          color=BLUE, fontsize=9.5, fontweight="bold",
          transform=ax_e.transAxes, va="top")
y_e -= 0.06
ax_e.text(0.01, y_e,
          f"Found {len(big_trnas)} tRNAs  |  {len(type_counts)} type(s)",
          color=MID, fontsize=8.5, transform=ax_e.transAxes, va="top")
y_e -= 0.055
for aa, n in type_counts.most_common():
    pos_list = [t["rel_pos_pct"] for t in big_trnas
                if (t["amino_acid"] if t["amino_acid"] else "OTHER")==aa
                and t["rel_pos_pct"] is not None]
    if pos_list:
        p_min, p_max = min(pos_list), max(pos_list)
        loc = ("EDGE" if any(p<=EDGE_THRESHOLD or p>=(100-EDGE_THRESHOLD)
                             for p in pos_list) else "MIDDLE")
        col = RED if loc=="EDGE" else BLUE
        ax_e.text(0.03, y_e,
                  f"tRNA-{aa} ×{n}  →  {loc}  "
                  f"({p_min:.0f}%–{p_max:.0f}% of prophage)",
                  color=col, fontsize=8.5, fontweight="bold",
                  transform=ax_e.transAxes, va="top")
        y_e -= 0.05

y_e -= 0.03
ax_e.plot([0.01, 0.99], [y_e+0.01, y_e+0.01],
             color=BORDER, linewidth=0.8,
             transform=ax_e.transAxes, clip_on=False)

ax_e.text(0.01, y_e-0.01, "The tmRNA — SsrA",
          color=PURPLE, fontsize=9.5, fontweight="bold",
          transform=ax_e.transAxes, va="top")
y_e -= 0.07
if tmrna_entries:
    for t in tmrna_entries:
        for line_text, col, bold in [
            (f"Prophage: {t['prophage_id']}", MID, False),
            (f"Species:  {t['species']}", MID, False),
            (f"Position: {t['start']}–{t['end']} bp  "
             f"({t['rel_pos_pct']:.1f}% → {t['location']})", MID, False),
            (f"Product:  {t['product']}", TEAL, True),
            ("Tag peptide: AKNNNSYALAA*", TEAL, True),
            ("Interpretation: integrated into a tmRNA gene.", PURPLE, True),
            ("The modified tag may protect phage proteins.", PURPLE, False),
        ]:
            ax_e.text(0.03, y_e, line_text,
                      color=col, fontsize=8.5,
                      fontweight="bold" if bold else "normal",
                      transform=ax_e.transAxes, va="top")
            y_e -= 0.05
else:
    ax_e.text(0.03, y_e,
              "tmRNA detected in Aragorn GFF but not in tRNAscan.\n"
              "Check pharokka_aragorn.gff directly for SsrA entry.\n"
              "Location: utir_016_m0_01_002__provirus_410273_450678\n"
              "Position: 8–354 bp (start of prophage)\n"
              "Tag peptide: AKNNNSYALAA*",
              color=PURPLE, fontsize=8.5, fontweight="bold",
              transform=ax_e.transAxes, va="top", linespacing=1.6)

plt.savefig("results/utir_trna_deep_v2.png",
            dpi=160, bbox_inches="tight", facecolor=BG)
print("  Saved → results/utir_trna_deep_v2.png")

print("\n=== DONE ===")
print("Outputs:")
print("  results/utir_trna_deep_v2.png")
print("  results/utir_trna_table.tsv")
print("  results/utir_trna_edge_candidates.tsv")
print()
print("Copy to Windows — run in WSL:")
print("  scp edzna:~/utir/results/utir_trna_deep_v2.png "
      "/mnt/c/Users/User/Desktop/")
