"""Summarise the within-lineage permutation null.

Only permutation outputs that pass a validity check are used: a run killed
part-way through leaves a large but incomplete file, which must not be counted
as a permutation.
"""
import numpy as np, pandas as pd, glob, sys
THR = 0.05 / 15618          # Bonferroni over unique presence-absence patterns
MIN_ROWS = 15000            # a complete pyseer run reports 16,151 variants
drugs = ["CIP", "CTX", "CTZ", "GEN"]

def valid(path):
    try:
        r = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    if "lrt-pvalue" not in r.columns or len(r) < MIN_ROWS or r["lrt-pvalue"].isna().any():
        return None
    return r

rows, strongest, dropped = [], [], 0
for dr in drugs:
    real = pd.read_csv(f"pyseer/lmm_{dr}.tsv", sep="\t")
    obs_n = int((real["lrt-pvalue"] < THR).sum())
    obs_min = real["lrt-pvalue"].min()
    counts, minima = [], []
    for f in sorted(glob.glob(f"perm/out/{dr}_p*.tsv")):
        r = valid(f)
        if r is None:
            dropped += 1
            continue
        counts.append(int((r["lrt-pvalue"] < THR).sum()))
        minima.append(r["lrt-pvalue"].min())
    counts, minima = np.array(counts), np.array(minima)
    rows.append(dict(drug=dr, observed=obs_n, perm_n=len(counts),
                     null_median=int(np.median(counts)), null_mean=round(counts.mean(), 1),
                     null_min=int(counts.min()), null_max=int(counts.max()),
                     null_95th=int(np.percentile(counts, 95)),
                     empirical_p=round((1 + (counts >= obs_n).sum()) / (1 + len(counts)), 3)))
    strongest.append(dict(drug=dr, perm_n=len(minima), observed_min_p=obs_min,
                          null_median_min_p=np.median(minima), null_best_min_p=minima.min(),
                          empirical_p=round((1 + (minima <= obs_min).sum()) / (1 + len(minima)), 3)))

S = pd.DataFrame(rows); S.to_csv("perm_summary.csv", index=False)
T = pd.DataFrame(strongest); T.to_csv("perm_strongest.csv", index=False)
if dropped:
    print(f"WARNING: {dropped} incomplete permutation files skipped", file=sys.stderr)
print("Number of significant loci per run:"); print(S.to_string(index=False))
print("\nStrongest single association per run:"); print(T.to_string(index=False))
