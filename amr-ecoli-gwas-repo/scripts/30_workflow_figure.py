"""Figure 1: study workflow (reviewer point 1)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({"font.family":"DejaVu Sans"})
LH=2.75; HEAD=7.0
C={"data":"#DCE9F5","qc":"#E8F1E4","struct":"#FFF2D9","pred":"#F6E3EC","assoc":"#EDE6F6","val":"#F1F1F1","out":"#FFFFFF"}
E={"data":"#3B6EA5","qc":"#4E8A3E","struct":"#9A6F1E","pred":"#A54A73","assoc":"#6A4FA0","val":"#555555","out":"#222222"}
def h(lines): return HEAD+LH*len(lines)
def box(ax,x,ytop,w,title,lines,k):
    H=h(lines); y=ytop-H
    ax.add_patch(FancyBboxPatch((x,y),w,H,boxstyle="round,pad=0.3,rounding_size=1.4",fc=C[k],ec=E[k],lw=1.2))
    ax.text(x+w/2,ytop-1.6,title,ha="center",va="top",fontsize=9,fontweight="bold",color=E[k])
    for i,l in enumerate(lines):
        ax.text(x+w/2,ytop-HEAD+0.4-i*LH,l,ha="center",va="top",fontsize=7.4)
    return y
def arrow(ax,x1,y1,x2,y2):
    ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",lw=1.1,color="#444",shrinkA=1,shrinkB=1))

B1=["Giessen dataset (Ren et al., 2022): public SNP matrix and phenotypes",
    "809 isolates with complete data for CIP, CTX, CTZ and GEN",
    "Variants called against E. coli K-12 MG1655 (≥20% call filter)"]
B2=["60,936 columns → 16,446 unique positions (exact duplicates collapsed)",
    "→ 16,435 polymorphic positions",
    "Stored codes N=0, A=1, G=2, C=3, T=4 (N = no variant called)",
    "Binary recoding: most frequent state = 0; all other states, incl. N = 1",
    "Sensitivity analysis: one-hot and original label encodings"]
B3=["Principal components (2–50) and genome-wide kinship matrix",
    "Maximum-likelihood phylogeny (FastTree) → patristic distances",
    "Lineages: 265 / 182 / 114 clusters at thresholds 0.05 / 0.1 / 0.2"]
B4a=["Resistant vs susceptible, per antibiotic",
     "Logistic regression, linear and RBF SVM,",
     "random forest",
     "Features: SNPs only, PCs only, SNPs + PCs",
     "Five-fold CV: random vs lineage-blocked",
     "(whole lineages held out)",
     "Out-of-fold ROC AUC"]
B4b=["Naive χ² (Bonferroni, BH-FDR)",
     "Linear model + 2–50 PCs",
     "Genomic inflation factor λ",
     "Linear mixed model with kinship (pyseer)",
     "Within-lineage CMH test",
     "Homoplasy (Fitch parsimony)",
     "Logistic refit: odds ratios, 95% CI"]
B5=["Within-lineage permutation null (200 per antibiotic)",
    "Synthetic-locus spike-in (odds ratios 2–16)",
    "Lineage-threshold sensitivity and clonality diagnostic",
    "Annotation against MG1655 (codon, synonymous / non-synonymous)",
    "Known resistance genes as positive controls; candidate loci excluding them"]
B6=["Lineage-level vs locus-level signal for each antibiotic"]
gap=4.5
tot=h(B1)+h(B2)+h(B3)+max(h(B4a),h(B4b))+h(B5)+h(B6)+5*gap+2
fig,ax=plt.subplots(figsize=(7.2,7.2*tot/100)); ax.set_xlim(0,100); ax.set_ylim(0,tot); ax.axis("off")
y=tot-1
y1=box(ax,7,y,86,"1  Data",B1,"data"); arrow(ax,50,y1,50,y1-gap+0.5)
y2=box(ax,7,y1-gap,86,"2  Preprocessing and genotype encoding",B2,"qc"); arrow(ax,50,y2,50,y2-gap+0.5)
y3=box(ax,7,y2-gap,86,"3  Population structure",B3,"struct")
arrow(ax,30,y3,26,y3-gap+0.5); arrow(ax,70,y3,74,y3-gap+0.5)
ya=box(ax,2,y3-gap,47,"4a  Prediction",B4a,"pred"); yb=box(ax,51,y3-gap,47,"4b  Association testing",B4b,"assoc")
y4=min(ya,yb); arrow(ax,25.5,ya,40,y4-gap+0.5); arrow(ax,74.5,yb,60,y4-gap+0.5)
y5=box(ax,7,y4-gap,86,"5  Validation and robustness",B5,"val"); arrow(ax,50,y5,50,y5-gap+0.5)
box(ax,17,y5-gap,66,"Outcome",B6,"out")
fig.savefig("Fig1_workflow.png",dpi=300,bbox_inches="tight",facecolor="white"); print("saved",tot)
