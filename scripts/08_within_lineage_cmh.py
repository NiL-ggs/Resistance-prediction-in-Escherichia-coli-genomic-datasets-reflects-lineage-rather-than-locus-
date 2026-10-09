import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy import stats
d=np.load("data.npz",allow_pickle=True); Gb,G,pos,Y=d["Gbin"].astype(float),d["G"],d["pos"],d["Y"]
D=np.load("tipdist.npy"); cl=fcluster(linkage(squareform(D,checks=False),"average"),0.1,"distance")
np.save("lineage_0.1.npy",cl)
drugs=["CIP","CTX","CTZ","GEN"]; m=Gb.shape[1]
def cmh(y):
    A=np.zeros(m);E=np.zeros(m);V=np.zeros(m)
    for c in np.unique(cl):
        s=cl==c; n=s.sum()
        if n<2: continue
        g=Gb[s]; yy=y[s]; r1=yy.sum(); c1=g.sum(0)
        a=yy@g; A+=a; E+=r1*c1/n; V+=r1*(n-r1)*c1*(n-c1)/(n*n*(n-1))
    with np.errstate(divide="ignore",invalid="ignore"): chi=(np.abs(A-E)-0.5)**2/V
    chi=np.nan_to_num(chi); p=stats.chi2.sf(chi,1); p[V<=0]=1; return p
def bh(p):
    o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1]; r=np.empty_like(q); r[o]=np.clip(q,0,1); return r
hom=pd.read_csv("homoplasy.csv").set_index("position")
out={}
for i,dr in enumerate(drugs):
    p=cmh(Y[:,i].astype(float)); q=bh(p); out[dr]=p
    o=np.argsort(p); r=int(np.where(pos[o]==3538947)[0][0])+1
    print(f"{dr}: CMH FDR<0.05: {(q<0.05).sum()}, Bonf: {(p<0.05/m).sum()} | 3538947 rank {r} p={p[pos==3538947][0]:.2e} | top: {list(pos[o[:5]])}")
pd.DataFrame({"position":pos,**{f"p_cmh_{k}":v for k,v in out.items()}}).to_csv("cmh_lineage.csv",index=False)
# lineage spread of top loci
for p in [3538947,4469454,4443644,4428463,4477553]:
    j=np.where(pos==p)[0][0]; v=Gb[:,j]==1
    lin=pd.Series(cl[v]).value_counts()
    print(p,"variant carriers:",v.sum(),"in",lin.size,"lineages; top lineage share",round(lin.iloc[0]/v.sum(),2),
          "| min changes on tree:",int(hom.loc[p,"min_changes_binary"]),"(genome median 13; pct",round((hom.min_changes_binary<hom.loc[p,'min_changes_binary']).mean()*100),")")
