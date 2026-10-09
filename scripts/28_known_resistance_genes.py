"""Reviewer point 5(i): behaviour of every matrix position lying in a curated set of
chromosomal genes / regulatory regions with established roles in resistance to the
four antibiotics.

Gene coordinates: from the GenBank record if present, otherwise the reconstructed table
(see scripts/annot_mg1655.py and reference/README_coordinates.md).

Outputs: known_gene_positions.csv, known_gene_summary.csv, known_gene_enrichment.csv
"""
import json, numpy as np, pandas as pd
from scipy import stats
from Bio.Seq import Seq

d=np.load("data.npz",allow_pickle=True); G,Gb,pos,maj,Y=d["G"],d["Gbin"],d["pos"],d["major"],d["Y"]
drugs=["CIP","CTX","CTZ","GEN"]; CODE={0:"N",1:"A",2:"G",3:"C",4:"T"}
import sys; sys.path.insert(0,"scripts"); import annot_mg1655 as A
genes=A.GENES; ref=A.REF; print("gene coordinates from:",A.SOURCE)
REGIONS={"marO (marR promoter/operator)":(1618909,1619119),"soxR–soxS intergenic":(4277384,4277468),
         "acrR–acrA intergenic":(485620,485760),"ampC promoter/attenuator":(4378945,4379006)}
KNOWN=[ # gene/region, mechanism class, antibiotics for which it is relevant, canonical sites
 ("gyrA","Quinolone target (QRDR)",["CIP"],"codons 83, 87"),("gyrB","Quinolone target",["CIP"],"codons 426, 447"),
 ("parC","Quinolone target (QRDR)",["CIP"],"codons 80, 84"),("parE","Quinolone target",["CIP"],"codons 416, 445, 458"),
 ("marR","Efflux regulation (mar)",["CIP","CTX","CTZ"],"inactivating"),("marA","Efflux regulation (mar)",["CIP","CTX","CTZ"],"expression"),
 ("marO (marR promoter/operator)","Efflux regulation (mar)",["CIP","CTX","CTZ"],"operator"),
 ("soxR","Efflux regulation (sox)",["CIP","CTX","CTZ"],"constitutive alleles"),("soxS","Efflux regulation (sox)",["CIP","CTX","CTZ"],"expression"),
 ("soxR–soxS intergenic","Efflux regulation (sox)",["CIP","CTX","CTZ"],"promoter"),
 ("rob","Efflux regulation",["CIP","CTX","CTZ"],"expression"),("acrR","Efflux regulation (acr)",["CIP","CTX","CTZ"],"inactivating"),
 ("acrR–acrA intergenic","Efflux regulation (acr)",["CIP","CTX","CTZ"],"promoter"),
 ("acrA","AcrAB–TolC efflux",["CIP","CTX","CTZ"],"—"),("acrB","AcrAB–TolC efflux",["CIP","CTX","CTZ"],"substrate-binding"),
 ("tolC","AcrAB–TolC efflux",["CIP","CTX","CTZ"],"—"),
 ("ompF","Porin",["CIP","CTX","CTZ"],"loss/altered"),("ompC","Porin",["CIP","CTX","CTZ"],"loss/altered"),
 ("ompR","Porin regulation",["CTX","CTZ"],"—"),("envZ","Porin regulation",["CTX","CTZ"],"—"),
 ("ampC","Chromosomal AmpC β-lactamase",["CTX","CTZ"],"—"),("ampC promoter/attenuator","AmpC hyperproduction",["CTX","CTZ"],"−42, −32"),
 ("ftsI","β-lactam target (PBP3)",["CTX","CTZ"],"YRIN/YRIK insertions"),
 ("fusA","Aminoglycoside (EF-G)",["GEN"],"—"),("cpxA","Aminoglycoside uptake (Cpx)",["GEN"],"—"),
 ("cpxR","Aminoglycoside uptake (Cpx)",["GEN"],"—"),("trkH","Aminoglycoside uptake",["GEN"],"—"),
]
lmm={}
for dr in drugs:
    t=pd.read_csv(f"pyseer/lmm_{dr}.tsv",sep="\t"); lmm[dr]=dict(zip(t.variant.str[4:].astype(int),t["lrt-pvalue"]))
assoc={dr:pd.read_csv(f"results/tables/{dr}_association_results_v2.csv").set_index("position") for dr in drugs}
cmh=pd.read_csv("cmh_lineage.csv").set_index("position")
hom=pd.read_csv("homoplasy.csv").set_index("position")
pi={p:i for i,p in enumerate(pos)}
THR_LMM=0.05/15618; THR_POS=0.05/len(pos)

def effect(p,gname):
    g=genes[gname]; j=pi[p]
    if g["strand"]=="+": off=p-g["start"]
    else: off=g["end"]-p
    codon_no=off//3+1; k=off%3
    cstart=g["start"]+3*(codon_no-1) if g["strand"]=="+" else g["end"]-3*codon_no+1
    cod=ref[cstart-1:cstart+2]
    cod=cod if g["strand"]=="+" else str(Seq(cod).reverse_complement())
    aa_ref=str(Seq(cod).translate())
    out=[]
    for code in sorted(set(G[:,j])-{0}):
        b=CODE[code] if g["strand"]=="+" else str(Seq(CODE[code]).complement())
        alt=cod[:k]+b+cod[k+1:]; aa=str(Seq(alt).translate())
        if alt==cod: continue
        out.append(f"{CODE[code]}:{aa_ref}{codon_no}{aa}"+(" (syn)" if aa==aa_ref else " (non-syn)"))
    return codon_no,"; ".join(out) if out else "reference state only"

rows=[]
for name,mech,rel,canon in KNOWN:
    if name in REGIONS: s,e=REGIONS[name]; ps=[p for p in pos if s<=p<=e]
    else: g=genes[name]; ps=[p for p in pos if g["start"]<=p<=g["end"]]
    if not ps:
        rows.append(dict(gene_region=name,mechanism=mech,relevant_to="/".join(rel),position=None)); continue
    for p in sorted(ps):
        j=pi[p]
        if name in REGIONS: cn,eff=None,"intergenic"
        else: cn,eff=effect(p,name)
        st=pd.Series(G[:,j]).map(CODE).value_counts()
        r=dict(gene_region=name,mechanism=mech,relevant_to="/".join(rel),position=int(p),codon=cn,effect=eff,
               states="; ".join(f"{a}:{n}" for a,n in st.items()),major=CODE[int(maj[j])],test_allele_freq=round(Gb[:,j].mean(),3),
               tree_changes=int(hom.loc[p,"min_changes_binary"]))
        for dr in rel:
            r[f"p_naive_{dr}"]=assoc[dr].loc[p,"p_naive"]; r[f"p_10PC_{dr}"]=assoc[dr].loc[p,"p_pc_adjusted"]
            r[f"p_LMM_{dr}"]=lmm[dr].get(p,np.nan); r[f"p_CMH_{dr}"]=cmh.loc[p,f"p_cmh_{dr}"]
        rows.append(r)
P=pd.DataFrame(rows); P.to_csv("known_gene_positions.csv",index=False)

summ=[]
for name,mech,rel,canon in KNOWN:
    sub=P[(P.gene_region==name)&P.position.notna()]
    r=dict(gene_region=name,mechanism=mech,relevant_to="/".join(rel),canonical_sites=canon,n_positions=len(sub),
           n_nonsyn=int(sub.effect.fillna("").str.contains("non-syn").sum()) if len(sub) else 0)
    for meth,thr in [("naive",THR_POS),("10PC",THR_POS),("LMM",THR_LMM),("CMH",THR_POS)]:
        ps=[sub[f"p_{meth}_{dr}"].min() for dr in rel] if len(sub) else []
        r[f"min_p_{meth}"]=np.nanmin(ps) if len(ps) else np.nan
        r[f"sig_{meth}"]=bool(np.nanmin(ps)<thr) if len(ps) else False
    summ.append(r)
S=pd.DataFrame(summ); S.to_csv("known_gene_summary.csv",index=False)
print(S[["gene_region","relevant_to","n_positions","n_nonsyn","min_p_naive","min_p_10PC","min_p_LMM","min_p_CMH","sig_LMM"]].to_string(index=False))

# enrichment: are positions in drug-relevant known genes ranked higher than the rest of the genome?
enr=[]
for dr in drugs:
    names=[n for n,_,rel,_ in KNOWN if dr in rel]
    kp=set(P[P.gene_region.isin(names)&P.position.notna()].position.astype(int))
    for meth in ["naive","10PC","LMM"]:
        if meth=="LMM": pv=pd.Series(lmm[dr])
        else: pv=assoc[dr]["p_naive" if meth=="naive" else "p_pc_adjusted"]
        inn=pv[pv.index.isin(kp)].dropna(); out=pv[~pv.index.isin(kp)].dropna()
        u=stats.mannwhitneyu(-np.log10(inn),-np.log10(out),alternative="greater")
        enr.append(dict(drug=dr,method=meth,n_known=len(inn),median_p_known=inn.median(),median_p_other=out.median(),
                        p_enrichment=u.pvalue))
E=pd.DataFrame(enr); E.to_csv("known_gene_enrichment.csv",index=False); print(E.to_string(index=False))
