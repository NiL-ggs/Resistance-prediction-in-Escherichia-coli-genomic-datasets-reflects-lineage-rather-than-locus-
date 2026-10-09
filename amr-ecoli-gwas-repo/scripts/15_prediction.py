import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
d=np.load("data.npz",allow_pickle=True); Gb,pos,Y=d["Gbin"].astype(float),d["pos"],d["Y"]
PCs=np.load("PCs50.npy"); lin=np.load("lineage_0.1.npy")
drugs=["CIP","CTX","CTZ","GEN"]; pi={p:i for i,p in enumerate(pos)}
thr=0.05/15618
rows=[]
for i,dr in enumerate(drugs):
    y=Y[:,i]
    t=pd.read_csv(f"pyseer/lmm_{dr}.tsv",sep="\t"); t["pos"]=t.variant.str[4:].astype(int)
    hits=[pi[p] for p in t.loc[t["lrt-pvalue"]<thr,"pos"] if p in pi]
    sets={
      "LMM loci only": Gb[:,hits] if hits else None,
      "Lineage (10 PCs) only": PCs[:,:10],
      "Lineage (20 PCs) only": PCs[:,:20],
      "LMM loci + 10 PCs": np.column_stack([Gb[:,hits],PCs[:,:10]]) if hits else None,
      "All 16,435 SNPs (RF)": Gb,
    }
    for name,X in sets.items():
        if X is None: rows.append(dict(drug=dr,model=name,n_features=0,auc_random_cv=np.nan,auc_lineage_cv=np.nan)); continue
        clf=RandomForestClassifier(n_estimators=300,random_state=0,n_jobs=2) if "RF" in name else LogisticRegression(max_iter=3000,C=1.0)
        # random 5-fold
        p1=cross_val_predict(clf,X,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),method="predict_proba")[:,1]
        # lineage-blocked 5-fold: whole lineages held out
        p2=cross_val_predict(clf,X,y,cv=GroupKFold(5),groups=lin,method="predict_proba")[:,1]
        rows.append(dict(drug=dr,model=name,n_features=X.shape[1],
                         auc_random_cv=round(roc_auc_score(y,p1),3),
                         auc_lineage_cv=round(roc_auc_score(y,p2),3)))
R=pd.DataFrame(rows); R.to_csv("prediction_auc.csv",index=False)
print(R.pivot(index="model",columns="drug",values=["auc_random_cv","auc_lineage_cv"]).to_string())
