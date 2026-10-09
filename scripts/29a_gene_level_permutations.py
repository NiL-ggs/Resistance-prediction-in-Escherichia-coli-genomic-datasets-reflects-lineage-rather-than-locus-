"""Within-lineage permuted phenotypes for the gene-level null (reviewer point 5ii).
Uses the same generator, seed and loop order as 16_permutation_setup.py, so the
first NPERM permutations here are identical to the first NPERM used for Figure 7."""
import numpy as np, pandas as pd, os
d=np.load("data.npz",allow_pickle=True); Y,S=d["Y"],d["samples"]
lin=np.load("lineage_0.1.npy"); rng=np.random.default_rng(0)
os.makedirs("perm_gene",exist_ok=True); NPERM=20
P=pd.DataFrame(index=[str(s) for s in S]); P.index.name="samples"
for k in range(NPERM):
    for i,dr in enumerate(["CIP","CTX","CTZ","GEN"]):
        y=Y[:,i].copy()
        for c in np.unique(lin):
            m=lin==c
            if m.sum()>1: y[m]=rng.permutation(y[m])
        P[f"{dr}_p{k}"]=y
P.to_csv("perm_gene/pheno_perm.tsv",sep="\t"); print(P.shape)
