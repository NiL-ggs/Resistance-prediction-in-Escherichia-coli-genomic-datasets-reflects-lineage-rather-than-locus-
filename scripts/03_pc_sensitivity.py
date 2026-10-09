import numpy as np, pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
d=np.load("data.npz",allow_pickle=True); Gb,pos,Y=d["Gbin"].astype(float),d["pos"],d["Y"]
drugs=["CIP","CTX","CTZ","GEN"]; n,m=Gb.shape
Gc=Gb-Gb.mean(0)
pca=PCA(50,random_state=0).fit(Gc); PCs=pca.transform(Gc); np.save("PCs50.npy",PCs)
print("cum var 1,2,5,10,20,50:",np.round(np.cumsum(pca.explained_variance_ratio_)[[0,1,4,9,19,49]],3))
def test(y,k):
    X=np.column_stack([np.ones(n)]+([PCs[:,:k]] if k else []))
    H=X@np.linalg.pinv(X.T@X)@X.T
    yr=y-H@y; Gr=Gb-H@Gb; ss=(Gr*Gr).sum(0); ok=ss>1e-8
    b=np.zeros(m); b[ok]=(Gr[:,ok].T@yr)/ss[ok]
    sse=((yr[:,None]-b*Gr)**2).sum(0); df=n-X.shape[1]-1
    t=np.zeros(m); t[ok]=b[ok]/np.sqrt(sse[ok]/df/ss[ok])
    p=2*stats.t.sf(abs(t),df); p[~ok]=1; return p
def bh(p):
    o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty_like(q); r[o]=np.clip(q,0,1); return r
def lam(p): return np.median(stats.chi2.isf(p,1))/stats.chi2.ppf(0.5,1)
rows=[]; P={}
for i,dr in enumerate(drugs):
    y=Y[:,i].astype(float)
    for k in [0,2,5,10,20,50]:
        p=test(y,k); P[(dr,k)]=p; q=bh(p)
        o=np.argsort(p); rank=int(np.where(pos[o]==3538947)[0][0])+1
        rows.append(dict(drug=dr,PCs=k,lambda_GC=round(lam(p),2),fdr_sig=int((q<0.05).sum()),bonf_sig=int((p<0.05/m).sum()),
            rank_3538947=rank,p_3538947=p[pos==3538947][0],top_hit=int(pos[o[0]])))
df=pd.DataFrame(rows); print(df.to_string(index=False)); df.to_csv("pc_sensitivity.csv",index=False)
np.save("P_pcs.npy",{f"{a}_{b}":v for (a,b),v in P.items()},allow_pickle=True)
