#!/bin/bash
cd "$HOME/utir" || exit 1
export GTDBTK_DATA_PATH=/home/hackllab/databases/gtdb/release226
echo "=== GTDB-Tk start $(date) ==="
/home/alex/bin/micromamba run -n gtdbtk26 gtdbtk classify_wf \
  --genome_dir data --extension fna \
  --out_dir results/gtdbtk_all \
  --skani_sketch_dir db/skani_sketch \
  --cpus 16 --pplacer_cpus 1 \
  > logs/gtdbtk_all.log 2>&1
echo "=== GTDB-Tk exit $? at $(date) ==="
