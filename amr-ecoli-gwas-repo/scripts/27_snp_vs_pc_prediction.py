"""Reviewer point 3: SNP-only vs PC-only vs SNP+PC prediction, under random and
lineage-blocked five-fold cross-validation.

PCs are computed in two ways. "fold" PCs (models PC, SNP+PC) are computed INSIDE
each training fold (PCA fitted on training isolates only, test isolates
projected), so no information from held-out isolates enters the features.
"global" PCs are the components of the full 809-isolate matrix (PCs50.npy, as
used in step 15); they use no phenotype information, but held-out isolates'
genotypes shape the axes, which is optimistic under lineage-blocked evaluation. Out-of-fold predictions are pooled to one AUC per model/scheme, and the
incremental value of SNP-level information over population structure is
summarised as AUC(SNP+PC) - AUC(PC only), with a 95% bootstrap interval over
isolates (paired, 1,000 resamples of the same out-of-fold predictions).

Outputs: snp_vs_pc_auc.csv, snp_vs_pc_delta.csv
"""
import sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import roc_auc_score

d=np.load("data.npz",allow_pickle=True); X=d["Gbin"].astype(float); Y=d["Y"]
lin=np.load("lineage_0.1.npy"); drugs=["CIP","CTX","CTZ","GEN"]
import os
RUN_DRUGS=os.environ.get("DRUGS","CIP,CTX,CTZ,GEN").split(",")
K=int(sys.argv[1]) if len(sys.argv)>1 else 20
CLFS=sys.argv[2].split(",") if len(sys.argv)>2 else ["LR","RF"]

def make_clf(name):
    if name=="LR": return LogisticRegression(max_iter=5000,C=1.0)
    return RandomForestClassifier(n_estimators=300,random_state=0,n_jobs=1)

PCg=np.load("PCs50.npy")[:,:K]
import os
MODELS=["SNP","PC","SNP+PC","PC (global)","SNP+PC (global)"]
if os.environ.get("CORE_ONLY"): MODELS=["SNP","PC","SNP+PC"]
if os.environ.get("MODELS"): MODELS=os.environ["MODELS"].split(",")
def features(model,tr,te):
    if model=="SNP": Xtr,Xte=X[tr],X[te]
    elif model=="PC (global)": Xtr,Xte=PCg[tr],PCg[te]
    elif model=="SNP+PC (global)": Xtr,Xte=np.hstack([X[tr],PCg[tr]]),np.hstack([X[te],PCg[te]])
    else:
        pca=PCA(K,random_state=0).fit(X[tr]); Ptr,Pte=pca.transform(X[tr]),pca.transform(X[te])
        if model=="PC": Xtr,Xte=Ptr,Pte
        else: Xtr,Xte=np.hstack([X[tr],Ptr]),np.hstack([X[te],Pte])
    return Xtr,Xte

def oof(model,clfname,y,splits):
    p=np.zeros(len(y))
    for tr,te in splits:
        Xtr,Xte=features(model,tr,te)
        if clfname=="LR":
            sc=StandardScaler().fit(Xtr); Xtr,Xte=sc.transform(Xtr),sc.transform(Xte)
            if model.startswith("SNP+PC"):
                # PCs enter as near-unpenalized covariates: rescaling the K standardized PC
                # columns by sqrt(m/K) gives the PC block the same total variance as the
                # SNP block, so the shared L2 penalty does not swamp 20 columns among 16,435
                w=np.sqrt(X.shape[1]/K); Xtr[:,-K:]*=w; Xte[:,-K:]*=w
        c=make_clf(clfname).fit(Xtr,y[tr]); p[te]=c.predict_proba(Xte)[:,1]
    return p

rng=np.random.default_rng(0); B=1000
rows=[]; deltas=[]
for i,dr in enumerate(drugs):
    if dr not in RUN_DRUGS: continue
    y=Y[:,i].astype(int)
    schemes={"random":list(StratifiedKFold(5,shuffle=True,random_state=0).split(X,y)),
             "lineage-blocked":list(GroupKFold(5).split(X,y,groups=lin))}
    boot=[rng.integers(0,len(y),len(y)) for _ in range(B)]
    for clfname in CLFS:
        for sch,spl in schemes.items():
            P={m:oof(m,clfname,y,spl) for m in MODELS}
            for m,p in P.items():
                rows.append(dict(drug=dr,classifier=clfname,n_pcs=K,cv=sch,model=m,auc=round(roc_auc_score(y,p),3)))
            for a,b in [("SNP+PC","PC"),("SNP","PC"),("SNP+PC","SNP")]:
                if a not in P or b not in P: continue
                obs=roc_auc_score(y,P[a])-roc_auc_score(y,P[b])
                bs=[roc_auc_score(y[ix],P[a][ix])-roc_auc_score(y[ix],P[b][ix]) for ix in boot if 0<y[ix].sum()<len(ix)]
                lo,hi=np.percentile(bs,[2.5,97.5])
                deltas.append(dict(drug=dr,classifier=clfname,n_pcs=K,cv=sch,contrast=f"{a} minus {b}",
                                   delta_auc=round(obs,3),ci_low=round(lo,3),ci_high=round(hi,3)))
            print(dr,clfname,sch,{m:round(roc_auc_score(y,p),3) for m,p in P.items()},flush=True)
tag=f"k{K}_{'_'.join(CLFS)}"+os.environ.get("TAG","")
pd.DataFrame(rows).to_csv(f"snp_vs_pc_auc_{tag}.csv",index=False)
pd.DataFrame(deltas).to_csv(f"snp_vs_pc_delta_{tag}.csv",index=False)
