#!/bin/bash
cd "$HOME/utir" || exit 1
MM=/home/alex/bin/micromamba

echo "=== GTDB-Tk  start $(date) ==="
"$MM" run -n gtdbtk26 gtdbtk classify_wf \
  --genome_dir data --extension fna \
  --out_dir results/gtdbtk_all \
  --cpus 16 --pplacer_cpus 1 \
  > logs/gtdbtk_all.log 2>&1
echo "=== GTDB-Tk  exit $? at $(date) ==="

echo "=== Pharokka  start $(date) ==="
"$MM" run -n pharokka18 pharokka.py \
  -i results/all_prophages_v113.fna \
  -o results/pharokka_v113 \
  -d /home/hackllab/databases/pharokka \
  -t 16 -m -f \
  > logs/pharokka_v113.log 2>&1
echo "=== Pharokka  exit $? at $(date) ==="
