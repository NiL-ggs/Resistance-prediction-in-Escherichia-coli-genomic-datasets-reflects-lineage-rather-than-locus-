import numpy as np, skbio, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
d=np.load("data.npz",allow_pickle=True); G,pos,Y,S=d["G"],d["pos"],d["Y"],d["samples"]
idx={s.replace("_"," "):i for i,s in enumerate(S)}
t=skbio.TreeNode.read("tree.nwk"); t=t.root_at_midpoint()
tips=list(t.tips()); y={}; 
for k,tp in enumerate(tips): y[id(tp)]=k
x={id(t):0.0}
for n in t.preorder(include_self=False): x[id(n)]=x[id(n.parent)]+(n.length or 0)
for n in t.postorder():
    if not n.is_tip(): y[id(n)]=np.mean([y[id(c)] for c in n.children])
fig=plt.figure(figsize=(11,13)); ax=fig.add_axes([0.03,0.1,0.55,0.84])
for n in t.preorder(include_self=False):
    p=n.parent; ax.plot([x[id(p)],x[id(n)]],[y[id(n)]]*2,c="k",lw=.3)
for n in t.postorder():
    if not n.is_tip(): ys=[y[id(c)] for c in n.children]; ax.plot([x[id(n)]]*2,[min(ys),max(ys)],c="k",lw=.3)
ax.set_ylim(-1,len(tips)); ax.axis("off")
order=[idx[tp.name] for tp in tips]
ax2=fig.add_axes([0.59,0.1,0.28,0.84]); ax2.set_ylim(-1,len(tips))
cols=[]
for i,dr in enumerate(["CIP","CTX","CTZ","GEN"]): cols.append((dr,[("#c44e52" if v else "#eeeeee") for v in Y[order,i]]))
def allele_col(p,cmap):
    j=np.where(pos==p)[0][0]; return [cmap[G[k,j]] for k in order]
cols.append(("3,538,947\n(yhgF)",allele_col(3538947,{0:"#f2b134",1:"#7a4fa3",2:"#bbbbbb",3:"#000",4:"#000"})))
cols.append(("4,469,454\n(mgtA)",allele_col(4469454,{0:"#2a9d8f",3:"#bbbbbb",1:"#000",2:"#000",4:"#000"})))
for c,(lab,colors) in enumerate(cols):
    ax2.bar([c]*len(colors),1,bottom=np.arange(len(colors))-0.5,color=colors,width=0.9,linewidth=0)
ax2.set_xticks(range(len(cols)),[c[0] for c in cols],rotation=90,fontsize=8); ax2.set_yticks([]); 
for s in ax2.spines.values(): s.set_visible(False)
fig.legend(handles=[Patch(color="#c44e52",label="Resistant"),Patch(color="#eeeeee",label="Susceptible"),
    Patch(color="#f2b134",label="3,538,947: N (ref T)"),Patch(color="#7a4fa3",label="3,538,947: A"),Patch(color="#bbbbbb",label="3,538,947 G / 4,469,454 C"),
    Patch(color="#2a9d8f",label="4,469,454: N (ref T)")],loc="lower center",ncol=3,fontsize=8,frameon=False)
fig.suptitle("Maximum-likelihood SNP phylogeny of 809 isolates (FastTree, midpoint-rooted)",fontsize=11)
plt.savefig("out/Fig_tree_rings.png",dpi=300); print("ok")
