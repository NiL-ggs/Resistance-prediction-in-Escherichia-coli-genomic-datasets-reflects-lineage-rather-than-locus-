import numpy as np, pandas as pd, glob, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
drugs=["CIP","CTX","CTZ","GEN"]; names=dict(CIP="Ciprofloxacin",CTX="Cefotaxime",CTZ="Ceftazidime",GEN="Gentamicin")
thr=0.05/15618
# ---- Figure: classifier comparison (main Figure 1 in new order)
R=pd.read_csv("classifier_comparison.csv")
order=["Logistic regression","SVM (linear)","SVM (RBF)","Random forest"]
cols=["#3a6ea5","#7aa6d0","#e07b39","#5a9367"]
fig,ax=plt.subplots(1,4,figsize=(15,3.8),sharey=True)
x=np.arange(2); w=0.19
for a,dr in zip(ax,drugs):
    s=R[R.drug==dr].set_index("model")
    for j,m in enumerate(order):
        a.bar(x+(j-1.5)*w,[s.loc[m,"auc_random"],s.loc[m,"auc_lineage_blocked"]],w,color=cols[j],label=m if dr=="CIP" else None)
    a.axhline(0.5,ls="--",c="k",lw=.8); a.set_xticks(x,["Random\n5-fold CV","Lineage-blocked\n5-fold CV"])
    a.set_title(names[dr]); a.set_ylim(0,1.0)
ax[0].set_ylabel("AUC")
fig.legend(loc="lower center",ncol=4,fontsize=9,frameon=False,bbox_to_anchor=(0.5,-0.06))
plt.tight_layout(); plt.savefig("out/Fig_classifiers.png",dpi=300,bbox_inches="tight"); plt.close()
# ---- Figure: permutation null (50 perms)
fig,ax=plt.subplots(1,4,figsize=(15,3.6))
S=pd.read_csv("perm_summary.csv").set_index("drug")
for a,dr in zip(ax,drugs):
    n=[int((pd.read_csv(f,sep="\t")["lrt-pvalue"]<thr).sum()) for f in sorted(glob.glob(f"perm/out/{dr}_p*.tsv"))]
    a.hist(n,bins=range(min(n),max(n)+2),color="#b9c6d6",edgecolor="white")
    a.axvline(S.loc[dr,"observed"],color="#c44e52",lw=2)
    a.set_title(f"{names[dr]}  (observed = {S.loc[dr,'observed']}, p = {S.loc[dr,'empirical_p']:.2f})",fontsize=9.5)
    a.set_xlabel("Significant loci per run")
ax[0].set_ylabel(f"Permutations (n={len(n)})")
fig.suptitle("Mixed-model hits: observed (red) against within-lineage permutation null",fontsize=11)
plt.tight_layout(); plt.savefig("out/Fig_permutation_null.png",dpi=300); plt.close()
print("ok")
