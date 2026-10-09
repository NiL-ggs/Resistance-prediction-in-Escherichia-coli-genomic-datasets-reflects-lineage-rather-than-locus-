#!/usr/bin/env bash
# pyseer LMM on each permuted phenotype; keeps only variants with p < 1e-4 to save space
set -euo pipefail
cd perm_gene
for k in $(seq 0 19); do for d in CIP CTX CTZ GEN; do
  c="${d}_p${k}"; [ -s "lmm_${c}.tsv" ] && continue
  pyseer --lmm --phenotypes pheno_perm.tsv --phenotype-column "$c" --pres ../pyseer/var.Rtab \
    --similarity ../pyseer/kinship.tsv --min-af 0.01 --max-af 0.99 --cpu 1 2>/dev/null \
    | awk -F'\t' 'NR==1 || $4<1e-4' > "lmm_${c}.tmp" && mv "lmm_${c}.tmp" "lmm_${c}.tsv"
done; done
echo FIN > done
