"""Lineage-threshold sensitivity: redefine lineages at 0.05 / 0.1 / 0.2 patristic
distance and repeat (a) blocked cross-validation and (b) permuted phenotypes."""
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
d=np.load("data.npz",allow_pickle=True); X=d["Gbin"].astype(float); Y=d["Y"]; S=d["samples"]
D=np.load("tipdist.npy"); Z=linkage(squareform(D,checks=False),"average")
drugs=["CIP","CTX","CTZ","GEN"]; rng=np.random.default_rng(0)
rows=[]
for thr in [0.05,0.1,0.2]:
    cl=fcluster(Z,thr,"distance"); np.save(f"lineage_{thr}.npy",cl)
    n_l=cl.max(); print(f"threshold {thr}: {n_l} lineages, largest {np.bincount(cl)[1:].max()}",flush=True)
    for i,dr in enumerate(drugs):
        y=Y[:,i]
        rf=RandomForestClassifier(n_estimators=300,random_state=0,n_jobs=1)
        p1=cross_val_predict(rf,X,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),method="predict_proba")[:,1]
        p2=cross_val_predict(rf,X,y,cv=GroupKFold(5),groups=cl,method="predict_proba")[:,1]
        rows.append(dict(threshold=thr,n_lineages=int(n_l),drug=dr,
                         auc_random=round(roc_auc_score(y,p1),3),
                         auc_blocked=round(roc_auc_score(y,p2),3),
                         drop=round(roc_auc_score(y,p1)-roc_auc_score(y,p2),3)))
        print("  ",dr,rows[-1]["auc_random"],rows[-1]["auc_blocked"],flush=True)
pd.DataFrame(rows).to_csv("threshold_cv.csv",index=False)
