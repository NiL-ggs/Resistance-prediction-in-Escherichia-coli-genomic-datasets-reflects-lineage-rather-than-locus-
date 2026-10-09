import numpy as np, pandas as pd, statsmodels.api as sm, warnings; warnings.filterwarnings("ignore")
exec(open("scripts/annotate_lib.py").read())
d=np.load("data.npz",allow_pickle=True); G,Gb,pos,Y=d["G"],d["Gbin"].astype(float),d["pos"],d["Y"]
PCs=np.load("PCs50.npy"); P=np.load("P_pcs.npy",allow_pickle=True).item()
drugs=["CIP","CTX","CTZ","GEN"]; lin=np.load("lineage_0.1.npy")
hom=pd.read_csv("homoplasy.csv").set_index("position"); cmh=pd.read_csv("cmh_lineage.csv").set_index("position")
lmm={}
for dr in drugs:
    t=pd.read_csv(f"pyseer/lmm_{dr}.tsv",sep="\t"); t["pos"]=t.variant.str[4:].astype(int); lmm[dr]=t.set_index("pos")
pi={p:i for i,p in enumerate(pos)}
def row(dr,p):
    i=drugs.index(dr); j=pi[p]; g=G[:,j]; al=[L[k] for k in np.unique(g) if k>0]; a=annotate(int(p),al)
    y=Y[:,i]; X=sm.add_constant(np.column_stack([Gb[:,j],PCs[:,:10]]))
    try: f=sm.Logit(y,X).fit(disp=0,maxiter=200); OR=np.exp(f.params[1]); lo,hi=np.exp(f.conf_int()[1]); pl=f.pvalues[1]
    except Exception: OR=lo=hi=pl=np.nan
    rates="; ".join(f"{'N(ref '+a['ref']+')' if k==0 else L[k]}: n={(g==k).sum()}, R={y[g==k].mean():.0%}" for k in np.unique(g))
    v=Gb[:,j]==1
    return {"Drug":dr,"Position (MG1655)":int(p),"Gene / context":a["gene"] if a["gene"]!="intergenic" else a["product"],
        "Product":a["product"] if a["gene"]!="intergenic" else "","Effect":a["changes"] or "intergenic",
        "Alleles: n, % resistant":rates,"Test-allele frequency":round(v.mean(),3),
        "p naive (χ²)":P[f"{dr}_0"][j],"p 10 PCs (linear)":P[f"{dr}_10"][j],"p 20 PCs (linear)":P[f"{dr}_20"][j],
        "OR 10 PCs (logistic)":OR,"OR 95% CI":f"{lo:.1f}–{hi:.1f}" if OR==OR else "","p 10 PCs (logistic)":pl,
        "p LMM (pyseer)":lmm[dr]["lrt-pvalue"].get(p,np.nan),"p within-lineage (CMH)":cmh.loc[p,f"p_cmh_{dr}"],
        "Min. changes on tree":int(hom.loc[p,"min_changes_binary"]),"Lineages carrying test allele":int(pd.Series(lin[v]).nunique())}
rows=[]
for dr in drugs:
    top=list(lmm[dr].sort_values("lrt-pvalue").index[:20])
    if 3538947 not in top: top.append(3538947)
    rows+= [row(dr,p) for p in top]
T=pd.DataFrame(rows)
# positive controls
pc=[]
for name,a,b,dr in [("gyrA",2336793,2339420,"CIP"),("parC",3163715,3165973,"CIP"),("parE",3173504,3175396,"CIP"),("marR",1619120,1619554,"CIP"),("ampC + 300 bp upstream",4377811,4379244,"CTX")]:
    m=np.where((pos>=a)&(pos<=b))[0]
    if len(m)==0: pc.append({"Drug":dr,"Region":name,"Position (MG1655)":"none in matrix"}); continue
    for j in m: r=row(dr,pos[j]); r["Region"]=name; pc.append(r)
PC=pd.DataFrame(pc)
S=pd.DataFrame([
 dict(Method="Naive χ² (no correction)",CIP=13688,CTX=11230,CTZ=1788,GEN=9010,Note="BH-FDR<0.05, 16,435 SNPs"),
 dict(Method="Linear model + 10 PCs",CIP=4624,CTX=1768,CTZ=2369,GEN=1201,Note="BH-FDR<0.05 (author's v2)")])
ps=pd.read_csv("pc_sensitivity.csv")
for k in [20,50]:
    s=ps[ps.PCs==k].set_index("drug").fdr_sig; S.loc[len(S)]=dict(Method=f"Linear model + {k} PCs",**s.to_dict(),Note="BH-FDR<0.05")
lm={dr:int((lmm[dr]["lrt-pvalue"]<0.05/15618).sum()) for dr in drugs}
S.loc[len(S)]=dict(Method="LMM (pyseer, SNP kinship)",**lm,Note="Bonferroni over 15,618 unique patterns (p<3.2e-6)")
cm={dr:int((cmh[f"p_cmh_{dr}"]<0.05/16435).sum()) for dr in drugs}
S.loc[len(S)]=dict(Method="Within-lineage CMH (182 tree lineages)",**cm,Note="Bonferroni p<3.0e-6")
with pd.ExcelWriter("out/GWAS_supplementary_tables.xlsx",engine="openpyxl") as w:
    pd.DataFrame({"Sheet":["S_summary","S_top_hits","S_positive_controls","S_PC_sensitivity"],
      "Contents":["Significant-SNP counts per method and drug","Top 20 LMM hits per drug (+ position 3,538,947) with annotation, alleles, effect sizes and p-values from every method",
                  "All SNPs in the matrix within known resistance genes; note gyrA codons 83/87 and parC are absent from the public matrix",
                  "Genomic inflation (λ), significant counts and rank of 3,538,947 for 0–50 PCs"],
      "Notes":["N = no variant called in Ren et al.'s matrix (MG1655 reference allele or absent region)","Test allele = all non-major codes pooled (author's v2 coding)","",""]}).to_excel(w,sheet_name="README",index=False)
    S.to_excel(w,sheet_name="S_summary",index=False); T.to_excel(w,sheet_name="S_top_hits",index=False); PC.to_excel(w,sheet_name="S_positive_controls",index=False); ps.to_excel(w,sheet_name="S_PC_sensitivity",index=False)
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
wb=load_workbook("out/GWAS_supplementary_tables.xlsx")
for ws in wb:
    for r in ws.iter_rows():
        for c in r:
            c.font=Font(name="Arial",size=10,bold=c.row==1)
            if c.row==1: c.fill=PatternFill("solid",fgColor="DDE6F0"); c.alignment=Alignment(wrap_text=True,vertical="top")
            if isinstance(c.value,float) and c.column_letter and ws.cell(1,c.column).value and str(ws.cell(1,c.column).value).startswith("p"): c.number_format="0.00E+00"
    for col in ws.columns: ws.column_dimensions[col[0].column_letter].width=min(45,max(10,max(len(str(c.value or "")) for c in col[:40])*0.9))
    ws.freeze_panes="A2"
wb.save("out/GWAS_supplementary_tables.xlsx")
print(S.to_string(index=False)); print(T[T["Position (MG1655)"]==3538947][["Drug","OR 10 PCs (logistic)","p LMM (pyseer)","p within-lineage (CMH)"]])
print(T.groupby("Drug").head(3)[["Drug","Position (MG1655)","Gene / context","Effect","p LMM (pyseer)"]].to_string(index=False))
print("nonsyn among top hits:",T[~T.Effect.str.contains("syn")&(T.Effect!="intergenic")][["Drug","Position (MG1655)","Gene / context","Effect"]].to_string(index=False))

# ---- export formatted table data for the document builders (13_*.js, 14_*.js) ----
import json
def _fmt(v):
    if isinstance(v,float):
        if np.isnan(v): return ""
        if abs(v)<0.001 and v!=0: return f"{v:.2e}"
        return f"{v:.3g}"
    return "" if v is None else str(v)
def _pack(df): return dict(cols=list(df.columns), rows=[[_fmt(v) for v in r] for r in df.values])
ps2=ps[["drug","PCs","lambda_GC","fdr_sig","bonf_sig","rank_3538947","p_3538947","top_hit"]].copy()
ps2.columns=["Drug","PCs","λ","FDR<0.05","Bonferroni","Rank of 3,538,947","p at 3,538,947","Top-ranked position"]
t3=T[["Drug","Position (MG1655)","Gene / context","Effect","Test-allele frequency","p naive (χ²)",
      "p 10 PCs (linear)","p LMM (pyseer)","p within-lineage (CMH)","OR 10 PCs (logistic)",
      "Min. changes on tree","Lineages carrying test allele"]].copy()
t3.columns=["Drug","Position","Gene / context","Effect","Test-allele freq.","p naive","p 10 PCs",
            "p LMM","p CMH","OR (10 PCs)","Changes on tree","Lineages"]
t4=PC[["Region","Drug","Position (MG1655)","Effect","Alleles: n, % resistant","p naive (χ²)",
       "p 10 PCs (linear)","p LMM (pyseer)"]].copy()
t4.columns=["Gene region","Drug","Position","Effect","Alleles: n, % resistant","p naive","p 10 PCs","p LMM"]
json.dump({"S1":_pack(S),"S2":_pack(ps2),"S3":_pack(t3),"S4":_pack(t4)}, open("tables.json","w"))
print("wrote tables.json")
