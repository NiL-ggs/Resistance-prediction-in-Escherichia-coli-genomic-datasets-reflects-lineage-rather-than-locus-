import numpy as np, pandas as pd, warnings, sys; warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
d=np.load("data.npz",allow_pickle=True); X=d["Gbin"].astype(float); Y=d["Y"]
lin=np.load("lineage_0.1.npy"); drugs=["CIP","CTX","CTZ","GEN"]
models={
 "Logistic regression": lambda: make_pipeline(StandardScaler(with_mean=False), LogisticRegression(max_iter=3000)),
 "SVM (linear)": lambda: make_pipeline(StandardScaler(with_mean=False), SVC(kernel="linear",probability=True,random_state=0)),
 "SVM (RBF)": lambda: make_pipeline(StandardScaler(with_mean=False), SVC(kernel="rbf",probability=True,random_state=0)),
 "Random forest": lambda: RandomForestClassifier(n_estimators=300,random_state=0,n_jobs=1),
}
dr=sys.argv[1]; i=drugs.index(dr); y=Y[:,i]; rows=[]
for name,f in models.items():
    p1=cross_val_predict(f(),X,y,cv=StratifiedKFold(5,shuffle=True,random_state=0),method="predict_proba")[:,1]
    p2=cross_val_predict(f(),X,y,cv=GroupKFold(5),groups=lin,method="predict_proba")[:,1]
    a1,a2=roc_auc_score(y,p1),roc_auc_score(y,p2)
    rows.append(dict(drug=dr,model=name,auc_random=round(a1,3),auc_lineage_blocked=round(a2,3),drop=round(a1-a2,3)))
    print(dr,name,round(a1,3),round(a2,3),flush=True)
pd.DataFrame(rows).to_csv(f"clf_{dr}.csv",index=False)
