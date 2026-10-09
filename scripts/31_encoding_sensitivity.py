"""Reviewer point 2: sensitivity of prediction to SNP encoding.
  binary   : most frequent state = 0, any other state (incl. N) = 1   (used throughout)
  one-hot  : one indicator per observed state per position (N, A, G, C, T kept separate;
             the most frequent state dropped as reference)
  label    : the stored codes N=0, A=1, G=2, C=3, T=4 used as numbers, as distributed
             (random forest only; an ordinal code has no meaning for a linear model)
Random vs lineage-blocked five-fold CV, as in step 20. Output: encoding_sensitivity.csv"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
d=np.load("data.npz",allow_pickle=True); G=d["G"]; Gb=d["Gbin"].astype(float); maj=d["major"]; Y=d["Y"]
lin=np.load("lineage_0.1.npy"); drugs=["CIP","CTX","CTZ","GEN"]
import os
RUN_DRUGS=os.environ.get("DRUGS","CIP,CTX,CTZ,GEN").split(",")
oh=[]
for k in range(5):
    m=(G==k).any(0)&(maj!=k)
    oh.append((G[:,m]==k).astype(np.float32))
OH=np.hstack(oh); print("one-hot columns:",OH.shape[1],flush=True)
enc={"binary":Gb,"one-hot":OH,"label (0-4)":G.astype(float)}
LR=lambda: make_pipeline(StandardScaler(with_mean=False),LogisticRegression(max_iter=3000))
RF=lambda: RandomForestClassifier(n_estimators=300,random_state=0,n_jobs=1)
# binary-encoding results for both classifiers are those of step 20 (identical folds);
# they are merged in from classifier_comparison.csv rather than recomputed
runs=[("Logistic regression","one-hot"),("Random forest","one-hot"),("Random forest","label (0-4)")]
rows=[]
for i,dr in enumerate(drugs):
    if dr not in RUN_DRUGS: continue
    y=Y[:,i]
    for cname,e in runs:
        X=enc[e]; f=LR if cname.startswith("Log") else RF
        p1=cross_val_predict(f(),X,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),method="predict_proba")[:,1]
        p2=cross_val_predict(f(),X,y,cv=GroupKFold(5),groups=lin,method="predict_proba")[:,1]
        a1,a2=roc_auc_score(y,p1),roc_auc_score(y,p2)
        rows.append(dict(drug=dr,classifier=cname,encoding=e,n_features=X.shape[1],auc_random=round(a1,3),auc_lineage_blocked=round(a2,3),drop=round(a1-a2,3)))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(f"encoding_part_{os.environ.get('DRUGS','all').replace(',','_')}.csv",index=False)
if os.environ.get("DRUGS"): raise SystemExit
cc=pd.read_csv("results/tables/classifier_comparison.csv")
cc=cc[cc.model.isin(["Logistic regression","Random forest"])].rename(columns={"model":"classifier"})
cc["encoding"]="binary"; cc["n_features"]=Gb.shape[1]
out=pd.concat([cc[["drug","classifier","encoding","n_features","auc_random","auc_lineage_blocked","drop"]],pd.DataFrame(rows)])
out.sort_values(["drug","classifier","encoding"]).to_csv("encoding_sensitivity.csv",index=False)
