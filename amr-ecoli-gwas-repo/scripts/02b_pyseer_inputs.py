"""Write pyseer inputs (var.Rtab, pheno.tsv, kinship.tsv) from data.npz.
Kinship = standardized genotype matrix Z, K = Z Z' / m."""
import numpy as np, pandas as pd
d=np.load("data.npz",allow_pickle=True); Gb=d["Gbin"].astype(float); pos=d["pos"]; S=[str(s) for s in d["samples"]]; Y=d["Y"]
pd.DataFrame(Gb.T.astype(int),index=[f"pos_{p}" for p in pos],columns=S).rename_axis("Gene").to_csv("pyseer/var.Rtab",sep="\t")
pd.DataFrame(Y,index=S,columns=["CIP","CTX","CTZ","GEN"]).rename_axis("samples").to_csv("pyseer/pheno.tsv",sep="\t")
Z=(Gb-Gb.mean(0))/Gb.std(0); K=Z@Z.T/Z.shape[1]
pd.DataFrame(K,index=S,columns=S).to_csv("pyseer/kinship.tsv",sep="\t")
print("written", Gb.shape)
