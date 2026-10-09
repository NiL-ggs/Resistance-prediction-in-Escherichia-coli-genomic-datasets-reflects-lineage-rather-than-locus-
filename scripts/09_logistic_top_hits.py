import numpy as np, pandas as pd, statsmodels.api as sm, warnings
warnings.filterwarnings("ignore")
from scipy import stats
d=np.load("data.npz",allow_pickle=True); G,Gb,pos,Y=d["G"],d["Gbin"].astype(float),d["pos"],d["Y"]
PCs=np.load("PCs50.npy")[:,:10]; drugs=["CIP","CTX","CTZ","GEN"]; L="NAGCT"
rows=[]
for i,dr in enumerate(drugs):
    r=pd.read_csv(f"results_v2/{dr}_association_results_v2.csv").head(50)
    loci=list(r.position)+([3538947] if 3538947 not in set(r.position) else [])
    y=Y[:,i]
    for rank,p in enumerate(loci,1):
        j=np.where(pos==p)[0][0]; X=sm.add_constant(np.column_stack([Gb[:,j],PCs]))
        try:
            f=sm.Logit(y,X).fit(disp=0,maxiter=200); b,se,pv=f.params[1],f.bse[1],f.pvalues[1]
            ci=np.exp([b-1.96*se,b+1.96*se])
        except Exception as e: b=se=pv=np.nan; ci=[np.nan,np.nan]
        g=G[:,j]; al=" / ".join(f"{L[k]}:{(g==k).sum()}" for k in range(5) if (g==k).sum())
        rows.append(dict(drug=dr,linear_rank=rank if rank<=50 else ">50",position=p,alleles=al,variant_freq=round(Gb[:,j].mean(),3),
            OR_logit=np.exp(b),CI_low=ci[0],CI_high=ci[1],p_logit=pv,p_linear=float(pd.read_csv(f"results_v2/{dr}_association_results_v2.csv").set_index("position").loc[p,"p_pc_adjusted"])))
df=pd.DataFrame(rows); df.to_csv("logistic_top50.csv",index=False)
for dr in drugs:
    s=df[df.drug==dr].dropna(subset=["p_logit"])
    rho=stats.spearmanr(s.p_linear,s.p_logit).correlation
    print(dr,"spearman(linear,logit) on top50:",round(rho,3)," n logit p<3.0e-6:",int((s.p_logit<0.05/16435).sum()),"/",len(s))
print(df[df.position==3538947][["drug","linear_rank","OR_logit","CI_low","CI_high","p_logit","p_linear"]].to_string(index=False))
