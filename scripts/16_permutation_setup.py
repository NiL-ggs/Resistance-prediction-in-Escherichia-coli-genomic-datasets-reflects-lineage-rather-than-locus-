import numpy as np, pandas as pd, os
d=np.load("data.npz",allow_pickle=True); Y,S=d["Y"],d["samples"]
lin=np.load("lineage_0.1.npy"); rng=np.random.default_rng(0)
os.makedirs("perm",exist_ok=True)
NPERM=20
P=pd.DataFrame(index=S); P.index.name="samples"
for k in range(NPERM):
    for i,dr in enumerate(["CIP","CTX","CTZ","GEN"]):
        y=Y[:,i].copy()
        for c in np.unique(lin):          # permute WITHIN lineages: destroys locus effect, keeps lineage effect
            m=lin==c
            if m.sum()>1: y[m]=rng.permutation(y[m])
        P[f"{dr}_p{k}"]=y
P.to_csv("perm/pheno_perm.tsv",sep="\t")
print(P.shape, "columns", list(P.columns[:4]))
