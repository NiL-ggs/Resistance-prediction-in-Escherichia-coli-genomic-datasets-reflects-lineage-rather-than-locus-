#!/usr/bin/env bash
set -u
for k in $(seq 0 199); do
  for d in CIP CTX CTZ GEN; do
    o="../perm/out/${d}_p${k}.tsv"
    [ -s "$o" ] && continue
    pyseer --lmm --phenotypes ../perm/pheno_perm.tsv --phenotype-column "${d}_p${k}" \
      --pres var.Rtab --similarity kinship.tsv --min-af 0.01 --max-af 0.99 --cpu 2 \
      > "$o" 2> "../perm/out/${d}_p${k}.log"
  done
done
touch ../perm/DONE200

# Re-run any output that is incomplete (a run killed part-way leaves a large but
# truncated file). 18_permutation_summary.py applies the same validity check.
python3 - <<'PY'
import pandas as pd, glob, os
bad=[f for f in glob.glob("../perm/out/*_p*.tsv")
     if (lambda r: r is None or "lrt-pvalue" not in r.columns or len(r)<15000
                   or r["lrt-pvalue"].isna().any())(
            pd.read_csv(f,sep="\t") if os.path.getsize(f)>0 else None)]
for f in bad: os.remove(f)
print(f"removed {len(bad)} incomplete files; re-run this script to regenerate them")
PY
