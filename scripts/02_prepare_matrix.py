import pandas as pd, numpy as np, re
from collections import defaultdict
geno = pd.read_csv("giessen/Giessen_dataset/cip_ctx_ctz_gen_multi_data.csv").set_index("prename")
pheno = pd.read_csv("giessen/Giessen_dataset/cip_ctx_ctz_gen_pheno.csv").set_index("prename")
common = geno.index.intersection(pheno.index); geno=geno.loc[common]; pheno=pheno.loc[common]
seen={}; cols=[]; pos=[]
for c in geno.columns:
    p=int(re.match(r"X(\d+)",c).group(1))
    if p not in seen: seen[p]=c; cols.append(c); pos.append(p)
G=geno[cols].values.astype(np.int8); pos=np.array(pos)
counts=np.stack([(G==k).sum(0) for k in range(5)])
major=counts.argmax(0)
Gbin=(G!=major).astype(np.int8)
vc=Gbin.sum(0); m=(vc>0)&(vc<len(G))
np.savez_compressed("data.npz", G=G[:,m], Gbin=Gbin[:,m], pos=pos[m], major=major[m],
    samples=np.array(common), Y=pheno[["CIP","CTX","CTZ","GEN"]].values)
print(G[:,m].shape, "order sorted?", np.all(np.diff(pos[m])>0))
