"""Positive control: spike a synthetic locus with a known effect into the real
genotype matrix and confirm the mixed model recovers it above the null."""
import numpy as np, pandas as pd, subprocess, os, sys
d=np.load("data.npz",allow_pickle=True); Gb,pos,Y,S=d["Gbin"],d["pos"],d["Y"],d["samples"]
lin=np.load("lineage_0.1.npy"); rng=np.random.default_rng(1)
os.makedirs("spike",exist_ok=True)
i=0  # ciprofloxacin
y0=Y[:,i].astype(int)
# synthetic locus: present in a random 30% of isolates, spread across lineages
# (so it is NOT a clone marker), with effect size OR applied by flipping phenotype
rows=[]
for OR in [2,4,8,16]:
    g=(rng.random(len(S))<0.30).astype(int)
    # build phenotype: baseline prevalence from real data, carriers get raised odds
    p0=y0.mean()
    odds=p0/(1-p0); pc=(odds*OR)/(1+odds*OR)
    y=np.where(g==1, (rng.random(len(S))<pc).astype(int), (rng.random(len(S))<p0).astype(int))
    rows.append((OR,g,y))
# write genotype matrix with the spiked locus appended, and phenotypes
o=np.argsort(pos); G=Gb[:,o]; P=pos[o]
var=pd.DataFrame(G.T,index=[f"pos_{p}" for p in P],columns=S)
ph=pd.DataFrame(index=S); ph.index.name="samples"
for OR,g,y in rows:
    var.loc[f"spike_OR{OR}"]=g
    ph[f"CIP_OR{OR}"]=y
var.rename_axis("Gene").to_csv("spike/var.Rtab",sep="\t")
ph.to_csv("spike/pheno.tsv",sep="\t")
print("written", var.shape, ph.shape)
