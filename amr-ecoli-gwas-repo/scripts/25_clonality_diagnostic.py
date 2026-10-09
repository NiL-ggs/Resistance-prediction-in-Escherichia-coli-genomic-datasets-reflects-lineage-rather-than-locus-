"""Diagnostic: how much of the phenotype's clonal structure does permutation at
each lineage threshold preserve? Used to decide which thresholds give a valid null."""
import numpy as np, pandas as pd
d=np.load("data.npz",allow_pickle=True); Y=d["Y"]
fine=np.load("lineage_0.05.npy")      # finest partition, used as the yardstick
drugs=["CIP","CTX","CTZ","GEN"]
def eta2(y,part):
    tot=((y-y.mean())**2).sum()
    within=sum(((y[part==c]-y[part==c].mean())**2).sum() for c in np.unique(part))
    return 1-within/tot
rows=[]
for i,dr in enumerate(drugs):
    y=Y[:,i].astype(float); r=dict(drug=dr, observed=round(eta2(y,fine),3))
    for t in ["0.05","0.1","0.2"]:
        cl=np.load(f"lineage_{t}.npy"); rng=np.random.default_rng(3); v=[]
        for _ in range(20):
            yp=y.copy()
            for c in np.unique(cl):
                m=cl==c
                if m.sum()>1: yp[m]=rng.permutation(yp[m])
            v.append(eta2(yp,fine))
        r[f"permuted_{t}"]=round(float(np.mean(v)),3)
    rows.append(r)
T=pd.DataFrame(rows); T.to_csv("clonality_diagnostic.csv",index=False)
print(T.to_string(index=False))
