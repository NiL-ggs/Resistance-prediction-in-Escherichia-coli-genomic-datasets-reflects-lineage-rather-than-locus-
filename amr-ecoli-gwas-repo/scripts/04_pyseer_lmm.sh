#!/usr/bin/env bash
# Linear mixed-model association with pyseer, using a kinship matrix built from
# the same SNPs. Inputs written by 02_prepare_matrix.py are expected in pyseer/.
set -euo pipefail
cd pyseer
for d in CIP CTX CTZ GEN; do
  pyseer --lmm \
    --phenotypes pheno.tsv --phenotype-column "$d" \
    --pres var.Rtab --similarity kinship.tsv \
    --min-af 0.01 --max-af 0.99 --cpu 2 \
    --output-patterns "pat_${d}.txt" > "lmm_${d}.tsv" 2> "log_${d}.txt"
done
# Bonferroni threshold = 0.05 / (number of unique presence-absence patterns)
wc -l < pat_CIP.txt
