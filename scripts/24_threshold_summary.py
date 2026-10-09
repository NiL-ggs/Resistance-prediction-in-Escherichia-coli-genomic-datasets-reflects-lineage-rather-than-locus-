"""Summarise the permutation null at alternative lineage-clustering thresholds."""
import numpy as np, pandas as pd, glob
THR=0.05/15618; drugs=["CIP","CTX","CTZ","GEN"]
real={d:int((pd.read_csv(f"pyseer/lmm_{d}.tsv",sep="\t")["lrt-pvalue"]<THR).sum()) for d in drugs}
rows=[]
for t in ["0.05","0.2"]:
    for d in drugs:
        n=[]
        for f in sorted(glob.glob(f"perm_thr/{t}_{d}_p*.tsv")):
            r=pd.read_csv(f,sep="\t")
            if "lrt-pvalue" in r.columns and len(r)>=15000 and not r["lrt-pvalue"].isna().any():
                n.append(int((r["lrt-pvalue"]<THR).sum()))
        n=np.array(n); obs=real[d]
        rows.append(dict(threshold=float(t),drug=d,perm_n=len(n),observed=obs,
            null_median=int(np.median(n)),null_range=f"{n.min()}–{n.max()}",
            empirical_p=round((1+(n>=obs).sum())/(1+len(n)),3)))
# add the main 0.1 analysis for comparison
main=pd.read_csv("perm_summary.csv")
for _,r in main.iterrows():
    rows.append(dict(threshold=0.1,drug=r.drug,perm_n=r.perm_n,observed=r.observed,
        null_median=r.null_median,null_range=f"{r.null_min}–{r.null_max}",empirical_p=r.empirical_p))
T=pd.DataFrame(rows).sort_values(["drug","threshold"])
T.to_csv("threshold_perm.csv",index=False)
print(T.to_string(index=False))
