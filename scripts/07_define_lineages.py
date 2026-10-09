import numpy as np, pandas as pd, skbio
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy import stats
d=np.load("data.npz",allow_pickle=True); Gb,pos,Y,S=d["Gbin"].astype(float),d["pos"],d["Y"],d["samples"]
t=skbio.TreeNode.read("tree.nwk")
dm=t.tip_tip_distances(); names=[n.replace(" ","_") for n in dm.ids]
order=[names.index(s) for s in S]; D=dm.data[np.ix_(order,order)]
Z=linkage(squareform(D,checks=False),"average")
for thr in [0.02,0.05,0.1,0.2]:
    cl=fcluster(Z,thr,"distance"); sz=np.bincount(cl)[1:]
    print(thr,"clusters",cl.max(),"size>=10:",(sz>=10).sum(),"singletons",(sz==1).sum(),"largest",sz.max())
np.save("tipdist.npy",D)
