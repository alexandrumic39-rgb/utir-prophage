# geNomad run provenance

## Superseded (do not use for results)
- results/subset_markeronly/  geNomad 1.9.9 (bioconda), DB v1.4, --disable-nn-classification
- results/genomad_long/       same config, 32 longitudinal isolates
- results/genomad_all/        same config, incomplete - killed partway
  Known problem: marker-only classification called a 64 kb ICE in
  utir_014_m0_04 a provirus (79 genes, 0 viral hallmarks). Rejected by the
  config below. Bioconda 1.9.9 also could not run the neural network at all
  (Keras 3 incompatibility).

## Primary
- results/genomad_v113/       geNomad 1.13.0 (pip), DB v1.9, full aggregated
                              classification (marker + neural network)

## Note
geNomad skips steps whose output files exist, ignoring database or parameter
changes. Always use a fresh output directory or --restart.

## vOTU clustering threshold audit (2026-10-31)

Catalogue: results/genomad_v113 (geNomad 1.13.0, DB v1.9, full
aggregated classification). 361 prophages >= 5 kb across 127 isolates.

Three pairwise rules tested, all at ANI >= 95%:
  min(AF_ref, AF_query) >= 85%                      -> 242 vOTUs
  max(AF_ref, AF_query) >= 85%                      -> 202 vOTUs
  max(AF) >= 85% AND len_short/len_long >= 0.70     -> 228 vOTUs  [ADOPTED]

Why min() was rejected: geNomad calls prophage boundaries
inconsistently between assemblies of the same strain. In participant
011 lineage 14 the same element was called at 34,572 bp in u1k1_07 and
27,088 bp in u2k2_07 (skani: ANI 100.00, AF 74.71 / 95.35). min() put
them in separate vOTUs, producing a spurious prophage loss AND gain
between her two UTI episodes. Reciprocal blastn showed all four
prophages present in both genomes at 100% identity and full coverage.

Why max() alone was rejected: contig-edge-truncated calls act as
bridges. A truncated fragment is contained in several distinct
full-length phages, and transitive union-find welds them together.
10 vOTUs had a shortest member under half their longest, the worst
spanning 14,137-133,906 bp across 19 members. None involved the 32
longitudinal isolates, which carry no edge-truncated prophages.

Why 0.70: the one verified real merge sits at ratio 0.78; every
chained cluster sits at 0.48 or below.

Primary matrix: results/lineage_prophage_matrix_v113_guarded.tsv
Sensitivity: _v113.tsv (min), _v113_maxaf.tsv (max, unguarded)
Scripts: prophage_compare.py, _maxaf.py, _v3.py
