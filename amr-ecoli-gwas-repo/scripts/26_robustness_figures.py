import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
drugs=["CIP","CTX","CTZ","GEN"]
names=dict(CIP="Ciprofloxacin",CTX="Cefotaxime",CTZ="Ceftazidime",GEN="Gentamicin")
THR=0.05/15618

# ---------- Figure: spike control ----------
sp=[]
for OR in [2,4,8,16]:
    r=pd.read_csv(f"spike/lmm_OR{OR}.tsv",sep="\t")
    s=r[r.variant==f"spike_OR{OR}"]
    sp.append((OR,float(s["lrt-pvalue"].iloc[0])))
obs_best={d:pd.read_csv(f"pyseer/lmm_{d}.tsv",sep="\t")["lrt-pvalue"].min() for d in drugs}
fig,ax=plt.subplots(figsize=(6.4,4.4))
x=[o for o,_ in sp]; y=[-np.log10(p) for _,p in sp]
ax.plot(x,y,"o-",color="#3a6ea5",lw=2,ms=8,label="Synthetic locus (known effect)")
ax.axhline(-np.log10(THR),ls="--",c="#c44e52",lw=1.2)
ax.text(6.5,-np.log10(THR)+1.2,"Significance threshold",color="#c44e52",fontsize=8.5)
ax.axhline(-np.log10(obs_best["CIP"]),ls=":",c="#666",lw=1.2)
ax.text(6.5,-np.log10(obs_best["CIP"])+1.2,"Strongest association in the observed data",color="#444",fontsize=8.5)
ax.set_xscale("log",base=2); ax.set_xticks(x,[str(v) for v in x])
ax.set_xlabel("Simulated odds ratio of the inserted locus"); ax.set_ylabel("−log10 p (mixed model)")
ax.set_title("Recovery of a synthetic locus inserted into the observed matrix",fontsize=10.5)
for xi,yi in zip(x,y): ax.annotate("rank 1",(xi,yi),textcoords="offset points",xytext=(8,6),fontsize=8,color="#3a6ea5")
ax.legend(fontsize=8.5,loc="upper left"); plt.tight_layout()
plt.savefig("out/Fig_spike.png",dpi=300); plt.close()

# ---------- Figure: threshold sensitivity ----------
cv=pd.read_csv("threshold_cv.csv"); pm=pd.read_csv("threshold_perm.csv"); cl=pd.read_csv("clonality_diagnostic.csv").set_index("drug")
fig,ax=plt.subplots(1,3,figsize=(14,4))
cols=dict(CIP="#3a6ea5",CTX="#e07b39",CTZ="#5a9367",GEN="#c44e52")
for d in drugs:
    s=cv[cv.drug==d].sort_values("threshold")
    ax[0].plot(s.threshold,s["drop"],"o-",color=cols[d],label=names[d])
ax[0].set_xlabel("Lineage-clustering threshold"); ax[0].set_ylabel("AUC lost under lineage-blocked CV")
ax[0].set_title("(a) Accuracy loss is present at every threshold",fontsize=10); ax[0].legend(fontsize=8)
for d in drugs:
    s=pm[pm.drug==d].sort_values("threshold")
    ax[1].plot(s.threshold,s.null_median,"o--",color=cols[d],alpha=.6)
    ax[1].plot(s.threshold,s.observed,"s-",color=cols[d],label=names[d])
ax[1].set_xlabel("Lineage-clustering threshold"); ax[1].set_ylabel("Significant loci")
ax[1].set_title("(b) Observed (solid) vs permutation null (dashed)",fontsize=10); ax[1].legend(fontsize=8)
for d in drugs:
    obs=cl.loc[d,"observed"]; vals=[cl.loc[d,f"permuted_{t}"] for t in ["0.05","0.1","0.2"]]
    ax[2].plot([0.05,0.1,0.2],vals,"o-",color=cols[d],label=names[d])
    ax[2].axhline(obs,ls=":",c=cols[d],alpha=.5)
ax[2].set_xlabel("Lineage-clustering threshold"); ax[2].set_ylabel("Phenotype clonality after permutation")
ax[2].set_title("(c) Dotted lines: observed clonality",fontsize=10); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig("out/Fig_threshold.png",dpi=300); plt.close()
print("ok")
