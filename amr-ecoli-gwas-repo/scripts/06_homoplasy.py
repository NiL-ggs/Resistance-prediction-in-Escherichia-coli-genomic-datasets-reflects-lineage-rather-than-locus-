import numpy as np, pandas as pd, skbio, sys
d=np.load("data.npz",allow_pickle=True); G,Gb,pos,S=d["G"],d["Gbin"],d["pos"],d["samples"]
t=skbio.TreeNode.read("tree.nwk"); t=t.root_at_midpoint() if hasattr(t,"root_at_midpoint") else t
t=t.bifurcate() if hasattr(t,"bifurcate") else t
idx={s.replace("_"," "):i for i,s in enumerate(S)}; idx.update({s:i for i,s in enumerate(S)})
def fitch(codes):  # codes: n x m small ints
    bits=(1<<codes).astype(np.uint8); m=codes.shape[1]; changes=np.zeros(m,np.int32); st={}
    for node in t.postorder():
        if node.is_tip(): st[id(node)]=bits[idx[node.name]]; continue
        ch=[st.pop(id(c)) for c in node.children]
        cur=ch[0]
        for c in ch[1:]:
            inter=cur&c; empty=inter==0; changes+=empty; cur=np.where(empty,cur|c,inter)
        st[id(node)]=cur
    return changes
ch_full=fitch(G.astype(np.int64)); ch_bin=fitch(Gb.astype(np.int64))
nst=np.array([len(np.unique(G[:,j])) for j in range(G.shape[1])])
out=pd.DataFrame(dict(position=pos,n_states=nst,min_changes_allele=ch_full,CI_allele=(nst-1)/ch_full,
                      min_changes_binary=ch_bin,CI_binary=1/ch_bin,variant_count=Gb.sum(0)))
out.to_csv("homoplasy.csv",index=False); print(out.describe().round(3).to_string())
print(len(list(t.tips())),"tips")
