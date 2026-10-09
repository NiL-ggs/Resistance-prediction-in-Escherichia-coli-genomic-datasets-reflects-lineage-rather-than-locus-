import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
d=np.load("data.npz",allow_pickle=True); pos,Y,Gb=d["pos"],d["Y"],d["Gbin"]
drugs=["CIP","CTX","CTZ","GEN"]; names=dict(CIP="Ciprofloxacin",CTX="Cefotaxime",CTZ="Ceftazidime",GEN="Gentamicin")
P=np.load("P_pcs.npy",allow_pickle=True).item(); PCs=np.load("PCs50.npy")
lmm={dr:pd.read_csv(f"pyseer/lmm_{dr}.tsv",sep="\t") for dr in drugs}
def lam(p): return np.median(stats.chi2.isf(np.clip(p,1e-300,1),1))/stats.chi2.ppf(.5,1)
# Fig QQ
fig,ax=plt.subplots(1,4,figsize=(16,4.2))
for a,dr in zip(ax,drugs):
    for lab,p,c in [("Naive",P[f"{dr}_0"],"#999"),("10 PCs",P[f"{dr}_10"],"#e07b39"),("20 PCs",P[f"{dr}_20"],"#c44e52"),("LMM (pyseer)",lmm[dr]["lrt-pvalue"].values,"#3a6ea5")]:
        p=np.sort(np.clip(p,1e-300,1)); e=-np.log10(np.arange(1,len(p)+1)/(len(p)+1))
        a.scatter(e,-np.log10(p),s=4,c=c,label=f"{lab} (λ={lam(p):.2f})")
    a.plot([0,5],[0,5],"k--",lw=.8); a.set_title(names[dr]); a.set_xlabel("Expected −log10 p"); a.set_ylim(0,80)
    a.legend(fontsize=7,markerscale=2)
ax[0].set_ylabel("Observed −log10 p"); plt.tight_layout(); plt.savefig("out/Fig_QQ.png",dpi=300); plt.close()
# PC sensitivity
ps=pd.read_csv("pc_sensitivity.csv")
fig,ax=plt.subplots(1,2,figsize=(10,4))
for dr in drugs:
    s=ps[ps.drug==dr]; ax[0].plot(s.PCs,s.lambda_GC,"o-",label=names[dr]); ax[1].plot(s.PCs,s.rank_3538947,"o-",label=names[dr])
ax[0].axhline(1,ls="--",c="k",lw=.8); ax[0].set_xlabel("Number of PCs"); ax[0].set_ylabel("Genomic inflation λ"); ax[0].set_yscale("log")
ax[1].set_yscale("log"); ax[1].invert_yaxis(); ax[1].set_xlabel("Number of PCs"); ax[1].set_ylabel("Rank of position 3,538,947 (log)")
ax[0].legend(); plt.tight_layout(); plt.savefig("out/Fig_PC_sensitivity.png",dpi=300); plt.close()
# PCA scatter
fig,ax=plt.subplots(1,4,figsize=(16,4))
for a,(i,dr) in zip(ax,enumerate(drugs)):
    a.scatter(PCs[:,0],PCs[:,1],c=np.where(Y[:,i]==1,"#c44e52","#8fb3d9"),s=8,alpha=.7)
    a.set_title(f"{names[dr]} (red = resistant)"); a.set_xlabel("PC1 (50.4%)")
ax[0].set_ylabel("PC2 (6.7%)"); plt.tight_layout(); plt.savefig("out/Fig_PCA.png",dpi=300); plt.close()
# phenotype correlation
fig,ax=plt.subplots(figsize=(4.2,3.6)); C=np.corrcoef(Y.T)
im=ax.imshow(C,cmap="Reds",vmin=0,vmax=1); ax.set_xticks(range(4),drugs); ax.set_yticks(range(4),drugs)
for i in range(4):
    for j in range(4): ax.text(j,i,f"{C[i,j]:.2f}",ha="center",va="center")
plt.colorbar(im,label="φ"); plt.title("Co-resistance"); plt.tight_layout(); plt.savefig("out/Fig_coresistance.png",dpi=300); plt.close()
print("pc2 var",np.round(np.var(PCs,0)[:2]/np.var(PCs,0).sum(),3))
