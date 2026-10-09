"""Add the revision sheets (Tables S8-S12) to the supplementary workbook built by step 10.

Reads the result CSVs from the working directory, or from results/tables/revision/ if they are
not there, and writes GWAS_supplementary_tables.xlsx in the working directory.
"""
import openpyxl, pandas as pd, numpy as np, json, shutil, os
from openpyxl.styles import Font
R="."; REV="results/tables/revision"
SRC="GWAS_supplementary_tables.xlsx"
if not os.path.exists(SRC): SRC=f"{R}/results/tables/GWAS_supplementary_tables.xlsx"
OUT="GWAS_supplementary_tables.xlsx"
if os.path.abspath(SRC)!=os.path.abspath(OUT): shutil.copy(SRC,OUT)
wb=openpyxl.load_workbook(OUT)
def find(name):
    for p in (f"{R}/{name}",f"{REV}/{name}"):
        if os.path.exists(p): return p
    return None
def add(name,df,width=18):
    if name in wb.sheetnames: del wb[name]
    ws=wb.create_sheet(name)
    ws.append(list(df.columns))
    for r in df.itertuples(index=False): ws.append([None if (isinstance(v,float) and np.isnan(v)) else v for v in r])
    for c in ws[1]: c.font=Font(name="Arial",size=10,bold=True)
    for row in ws.iter_rows(min_row=2):
        for c in row: c.font=Font(name="Arial",size=10)
    for i,col in enumerate(df.columns,1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width=width
    ws.freeze_panes="A2"
    return ws
# S8 — within-lineage CMH, per position
cmh=pd.read_csv(find("cmh_lineage.csv") or f"{R}/results/tables/cmh_lineage.csv")
add("S_cmh",cmh,16)
# S9 — prediction from SNPs / PCs / both
import glob
pre=f"{REV}/feature_set_auc.csv"
parts=sorted(glob.glob(f"{R}/snp_vs_pc_auc_*.csv"))
if parts:
    A=pd.concat([pd.read_csv(f) for f in parts],ignore_index=True).drop_duplicates(["drug","classifier","n_pcs","cv","model"])
    A["pca_fitted_on"]=np.where(A.model.str.contains("global"),"all 809 isolates","training fold")
    A["model"]=A.model.str.replace(" (global)","",regex=False)
    A=A[["drug","classifier","n_pcs","cv","model","pca_fitted_on","auc"]]
elif os.path.exists(pre): A=pd.read_csv(pre)
else: A=None
if A is not None: add("S_feature_sets",A,18)
dci=find("snp_vs_pc_delta_k20_LR_unweighted.csv") or (f"{REV}/feature_set_delta_ci.csv" if os.path.exists(f"{REV}/feature_set_delta_ci.csv") else None)
if dci: add("S_feature_sets_CI",pd.read_csv(dci),18)
# S10 — encoding
for sheet,fn in [("S_encoding","encoding_sensitivity.csv"),
                 ("S_known_genes","known_gene_summary.csv"),
                 ("S_known_gene_positions","known_gene_positions.csv"),
                 ("S_known_gene_enrichment","known_gene_enrichment.csv"),
                 ("S_candidate_genes","candidate_genes.csv"),
                 ("S_candidate_loci","candidate_loci.csv"),
                 ("S_gene_level_null","gene_level_null.csv")]:
    p=find(fn)
    if p: add(sheet,pd.read_csv(p),18)
# README rows
rd=wb["README"]
new=[("S_cmh","Table S8 — within-lineage Cochran–Mantel–Haenszel p-value at every position, all four antibiotics","Stratified by the 182 lineages at threshold 0.1"),
 ("S_feature_sets","Table S9 — AUC for SNPs only, principal components only and both, under random and lineage-blocked folds","Components computed within training folds unless stated"),
 ("S_feature_sets_CI","Table S9 — differences between feature sets with 95% paired bootstrap intervals","1,000 resamples of isolates"),
 ("S_encoding","Table S10 — classification accuracy under binary, one-hot and integer-code encodings","Same folds as Table S1"),
 ("S_known_genes","Table S11 — known resistance genes and regions: positions present and smallest p-value per method","27 genes and regulatory regions"),
 ("S_known_gene_positions","Table S11 — per-position detail for the known resistance genes","Allele states, effects and p-values"),
 ("S_known_gene_enrichment","Table S11 — enrichment of known-gene positions among strong associations","One-sided Mann–Whitney on −log10 p"),
 ("S_candidate_genes","Table S12 — candidate loci outside known resistance genes, grouped by gene or interval","Includes permutation recurrence"),
 ("S_candidate_loci","Table S12 — per-position detail for the candidate loci","Within-lineage p, homoplasy, lineage concentration"),
 ("S_gene_level_null","Table S12 — gene-level significance against the within-lineage permutation null","Šidák within gene, Bonferroni across genes")]
have={r[0].value for r in rd.iter_rows(min_row=2)}
for row in new:
    if row[0] in have: continue
    rd.append(list(row))
    for c in rd[rd.max_row]: c.font=Font(name="Arial",size=10)
# fix the two renumbered main-text references in README notes
for r in rd.iter_rows(min_row=2):
    if r[2].value and "main-text Figure 1" in r[2].value: r[2].value=r[2].value.replace("Figure 1","Figure 2")
    if r[2].value and "main-text Figure 7" in r[2].value: r[2].value=r[2].value.replace("Figure 7","Figure 6")
    if r[2].value and "main-text Table 2" in r[2].value: r[2].value=r[2].value.replace("Table 2","Table 3")
wb.save(OUT); print("wrote",OUT,"with",len(wb.sheetnames),"sheets")
