# UTIr Prophage Study — Key Findings

## Overview
127 bacterial genomes from 21 women (recurrent UTI cohort, UMC Utrecht).  
361 prophages detected by geNomad 1.13.0 across 47 bacterial species.  
All analyses run on Edzna (24 CPUs, 62 GB RAM).

---

## Finding 1 — Prophages are nearly universal in urinary bacteria
- **85.8%** of isolates (109/127) carry at least one prophage
- Mean **2.85 prophages** per isolate; maximum **16** (E. coli, participant 004)
- Prevalence is high across nearly every species — not driven by one outlier
- Matches the only comparable published value (86% in 181 urinary isolates)

## Finding 2 — Prophage content is fixed within a persisting strain
- 9 persisting lineages tracked across 4 women, 32 isolates total
- **18 of 19** tracked prophages completely unchanged between timepoints
- One satellite lost between two cultures of one UTI episode (032/6)
- Lineage 032/2 (L. jensenii): 8 colonies, 3 timepoints, all carry exactly
  3 prophages at 8.02–8.05% genomic fraction (range 0.03%)
- Prophage load does **not** differ between UTI-episode and UTI-free visits
  once host species is accounted for — clean negative result

## Finding 3 — Prophages are extremely diverse and mostly unknown
- 228 distinct virus types (vOTUs); **170 found only once**
- **88%** receive no genus classification from the INPHARED database
- Named genera match their expected hosts (Phifelvirus in enterococci etc.)
  — sanity check on the whole catalogue

## Finding 4 — The Sit iron transport system is the strongest virulence link
- Sit genes found in **7 isolates** across participants 004, 007, 011, 014, 025
- **Complete sitABCD operon** in 3 isolates (011_u1k1_04, 025_m0_05/06)
- Iron is severely restricted in urine — Sit gives bacteria a competitive
  advantage in the iron-poor urinary environment
- Complete operon already present at **healthy baseline** in participant 025
  — predates any UTI, suggesting it contributes to infection establishment

## Finding 5 — Complete CDT toxin operon in participant 014 (CRITICAL)
- cdtA (93.1%), cdtB (79.3% PHROGs), cdtC (88.8%) all in **one prophage**
  (`utir_014_u1k1_05_001__provirus_2549520_2581857`, CDS_0011/0012/0013)
- cdtB was not detected by VFDB (sparse reference entries) but is
  unambiguously annotated by Pharokka as "cytolethal distending toxin B subunit"
- CDT is an AB2 toxin: binding subunits (CdtA, CdtC) deliver the DNase
  (CdtB) into human immune cells, causing DNA double-strand breaks,
  cell cycle arrest, and immune cell death
- Participant 014 had **4 UTI episodes** — most in the cohort
- Also carries: emrE, baeR, partial Sit, nmpC across other isolates
- Most virulence-diverse participant in the dataset

## Finding 6 — Bile salt hydrolase in all E_B. lactis isolates
- BSH gene (82.5% identity) in all 4 E_B. lactis isolates from
  participants 011, 025, 029 — all copies of the **same shared prophage**
- BSH is conserved across all copies → under purifying selection → functional
- Likely a relic of gut ancestry maintained for stress survival in urine

## Finding 7 — Resistance genes are rare (effectively negative result)
- **4 CARD hits** in 17,129 prophage genes (0.023%)
- emrE (efflux pump) ×2, baeR (efflux regulator) ×2
- These prophages are NOT spreading antibiotic resistance

## Finding 8 — tRNA analysis reveals integration sites and translation tools
- **180 tRNAs + 1 tmRNA** across 110 prophages (30%)
- **63 edge-positioned** (integration site candidates) — Arg tRNA dominant
- **114 middle-positioned** — codon-bias correction tools
- S. vaginalis 81 kb phage carries **14 distinct tRNA types** in one central
  block — host-range expansion strategy (novel in urinary tract phages)
- tmRNA at position 8 bp in utir_016 prophage with modified tag peptide
  AKNNNSYALAA — extremely rare, integration into tmRNA gene
- tRNA-carrying prophages score High-quality on CheckV 67% vs 46%
  for those without — independent validation of edge-tRNA integration

## Finding 9 — terL phylogeny reveals Sit gene evolution
- Sit genes exist in **at least 2 completely unrelated phage lineages**
  (distance ~1.17 between groups) — independent acquisition events
- Shared prophages between women are **co-inherited with bacterial strains**
  not recently transferred (025/029 E_B. lactis: 0.000 distance;
  006/010 E. faecalis: 0.000 distance)

## Finding 10 — Prophage sharing follows strain sharing
- 8 of 19 multi-isolate strain clusters span more than one woman
- Shared prophages track shared bacterial strains exactly
- No evidence for direct phage transfer between unrelated bacteria

---

## Status
- **Descriptive analysis: complete**
- **Clinical comparison (cases vs controls): blocked on metadata**
  (which participants are recurrent-UTI cases vs controls not provided)

---

## Scripts
| Script | Purpose |
|--------|---------|
| ani_clusters.py | Species and strain ANI clustering |
| prophage_compare_v3.py | vOTU clustering, lineage stability |
| classify_elements.py | Prophage classification (phage/satellite/uncertain) |
| utir_virulence.py | VFDB virulence gene analysis |
| utir_trna_analysis_v4.py | tRNA integration site and codon bias analysis |
| utir_phylogeny_v3.py | terL phylogenetic tree |
| utir_followups.py | Load per Mb, sizes, terL, cargo, edge analysis |
| utir_followups2.py | Taxonomy, twins, large elements, reconciliation |
| utir_brief_gaps.py | UTI vs UTI-free, genomic fraction |
| utir_adhoc.py | Sampling structure, edge cross-tab, fragmentation |

## Key result files
| File | Contents |
|------|---------|
| results/utir_terl_phylogeny_v3.png | terL NJ tree coloured by species |
| results/utir_terl_distances_v3.tsv | Pairwise terL distances |
| results/utir_virulence.png | 6-panel virulence figure |
| results/utir_virulence_table.tsv | All 41 VFDB hits with metadata |
| results/utir_trna_deep_v2.png | tRNA position and type analysis |
| results/utir_trna_table.tsv | All 180 tRNAs with EDGE/MIDDLE classification |
| results/utir_trna_edge_candidates.tsv | 63 integration site candidates |
