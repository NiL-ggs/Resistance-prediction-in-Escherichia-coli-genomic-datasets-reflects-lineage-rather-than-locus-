"""Reviewer point 5(ii): candidate loci outside known resistance genes.
(a) Position level: every mixed-model-significant position (p < 0.05/15,618) not in a
    known resistance gene/region (step 28), grouped by gene, with within-lineage (CMH)
    support, homoplasy, lineage concentration of carriers, and the frequency with which
    the same gene carries a significant position under within-lineage permutation.
(b) Gene level: min-p per gene with Sidak correction for the number of positions in the
    gene, Bonferroni across genes; counts compared with the permutation null.
Outputs: candidate_loci.csv, candidate_genes.csv, gene_level_null.csv"""
import glob, os, numpy as np, pandas as pd
import sys; sys.path.insert(0,"scripts"); import annot_mg1655 as A
d=np.load("data.npz",allow_pickle=True); G,Gb,pos=d["G"],d["Gbin"],d["pos"]
lin=np.load("lineage_0.1.npy"); drugs=["CIP","CTX","CTZ","GEN"]; pi={p:i for i,p in enumerate(pos)}
THR=0.05/15618
known=set(pd.read_csv("known_gene_positions.csv").position.dropna().astype(int))
cmh=pd.read_csv("cmh_lineage.csv").set_index("position"); hom=pd.read_csv("homoplasy.csv").set_index("position")
gene_of={p:(A.gene_at(p) or A.flanks(p)) for p in pos}
# Leading positions were annotated from the GenBank record in step 10 (Table S4); use those
# labels for reporting where available, the coordinate helper otherwise.
import json
_tj="tables.json" if os.path.exists("tables.json") else "results/tables/tables.json"
S4={int(r[1]):(r[2],r[3]) for r in json.load(open(_tj))["S3"]["rows"]}
def lmm(path):
    t=pd.read_csv(path,sep="\t"); return dict(zip(t.variant.str[4:].astype(int),t["lrt-pvalue"]))
obs={dr:lmm(f"pyseer/lmm_{dr}.tsv") for dr in drugs}
perms={dr:[lmm(f) for f in sorted(glob.glob(f"perm_gene/lmm_{dr}_p*.tsv"))] for dr in drugs}
npm={dr:len(v) for dr,v in perms.items()}; print("permutations:",npm)
# genes carrying >=1 significant position in each permutation
perm_genes={dr:[{gene_of[p] for p,pv in r.items() if pv<THR and p not in known} for r in perms[dr]] for dr in drugs}
rows=[]
for dr in drugs:
    for p,pv in obs[dr].items():
        if pv>=THR or p in known: continue
        j=pi[p]; g,eff=A.annotate(p,G[:,j]); car=Gb[:,j]==1
        if p in S4: g,eff=S4[p]; eff=eff if eff!="intergenic" else "intergenic"
        else: g=gene_of[p]
        ls=pd.Series(lin[car]).value_counts()
        rows.append(dict(drug=dr,position=p,gene=g,effect=eff,test_allele_freq=round(car.mean(),3),
             p_LMM=pv,p_CMH=cmh.loc[p,f"p_cmh_{dr}"],tree_changes=int(hom.loc[p,"min_changes_binary"]),
             carrier_lineages=int(ls.size),top_lineage_share=round(ls.iloc[0]/car.sum(),2),
             in_4_42_4_60Mb=bool(4_420_000<=p<=4_600_000),
             gene_hit_freq_in_permutations=round(np.mean([gene_of[p] in s for s in perm_genes[dr]]),2)))
L=pd.DataFrame(rows).sort_values(["drug","p_LMM"]); L.to_csv("candidate_loci.csv",index=False)
def nonsyn(eff):
    import re
    return any(re.search(r':[A-Z*]\d+[A-Z*]',x) and "(syn)" not in x for x in str(eff).split("; "))
# gene summary
gs=[]
for g,sub in L.groupby("gene"):
    gs.append(dict(gene=g,drugs="/".join(sorted(sub.drug.unique(),key=drugs.index)),n_drugs=sub.drug.nunique(),
        positions=", ".join(str(x) for x in sorted(sub.position.unique())),
        any_non_syn=any(nonsyn(e) for e in sub.effect),best_p_LMM=sub.p_LMM.min(),
        within_lineage_support=bool((sub.p_CMH<0.05/len(pos)).any()),min_tree_changes=int(sub.tree_changes.min()),
        max_top_lineage_share=sub.top_lineage_share.max(),in_4_42_4_60Mb=bool(sub.in_4_42_4_60Mb.all()),
        max_perm_hit_freq=sub.gene_hit_freq_in_permutations.max()))
S=pd.DataFrame(gs).sort_values(["n_drugs","best_p_LMM"],ascending=[False,True]); S.to_csv("candidate_genes.csv",index=False)
pd.set_option("display.width",250); print(S.to_string(index=False))
# gene-level Sidak test
genes={}
for p in pos:
    if p in known: continue
    genes.setdefault(gene_of[p],[]).append(p)
def gene_sig(pv):
    out=[]
    for g,ps in genes.items():
        v=[pv[p] for p in ps if p in pv]
        if not v: continue
        k=len(ps); pg=1-(1-min(v))**k; out.append((g,pg))
    return out
nul=[]
for dr in drugs:
    o=gene_sig(obs[dr]); ng=len(o); thr=0.05/ng
    nobs=sum(pg<thr for _,pg in o)
    nperm=[sum(pg<thr for _,pg in gene_sig(r)) for r in perms[dr]]
    nul.append(dict(drug=dr,genes_tested=ng,threshold=thr,observed_sig_genes=nobs,null_median=np.median(nperm),
        null_min=min(nperm),null_max=max(nperm),empirical_p=(1+sum(x>=nobs for x in nperm))/(1+len(nperm)),
        sig_genes=", ".join(g for g,pg in sorted(o,key=lambda x:x[1]) if pg<thr)))
N=pd.DataFrame(nul); N.to_csv("gene_level_null.csv",index=False); print(N.to_string(index=False))
