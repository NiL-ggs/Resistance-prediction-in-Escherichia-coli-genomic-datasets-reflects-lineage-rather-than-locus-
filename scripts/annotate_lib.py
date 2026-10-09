from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np, bisect
r=SeqIO.read("reference/NC_000913.3.gb","genbank"); SEQ=str(r.seq)
cds=[]
for f in r.features:
    if f.type=="CDS":
        q=f.qualifiers; cds.append((int(f.location.start),int(f.location.end),f.location.strand,
            q.get("gene",q.get("locus_tag",["?"]))[0],q.get("product",[""])[0],f))
genes=[(int(f.location.start),int(f.location.end),f.location.strand,f.qualifiers.get("gene",f.qualifiers.get("locus_tag",["?"]))[0])
       for f in r.features if f.type=="gene"]
genes.sort(); starts=[g[0] for g in genes]
L="NAGCT"
def annotate(p, alts=()):
    ref=SEQ[p-1]; hits=[c for c in cds if c[0]<p<=c[1]]
    if hits:
        s,e,strand,g,prod,f=hits[0]
        if strand==1: off=p-1-s
        else: off=e-p
        cn=off//3+1; ci=off%3
        codon_start = s+ (cn-1)*3 if strand==1 else None
        if strand==1: cod=SEQ[s+(cn-1)*3:s+cn*3]
        else: cod=str(Seq(SEQ[e-cn*3:e-(cn-1)*3]).reverse_complement())
        ch=[]
        for a in alts:
            b=a if strand==1 else str(Seq(a).reverse_complement())
            nc=cod[:ci]+b+cod[ci+1:]
            aa0=str(Seq(cod).translate()); aa1=str(Seq(nc).translate())
            ch.append(f"{a}:{aa0}{cn}{aa1}" + (" (syn)" if aa0==aa1 else ""))
        return dict(ref=ref,gene=g,product=prod,context="CDS",codon=cn,changes="; ".join(ch))
    i=bisect.bisect(starts,p)
    left=genes[i-1] if i>0 else None; right=genes[i] if i<len(genes) else None
    desc=f"between {left[3]}({'+' if left[2]==1 else '-'}) and {right[3]}({'+' if right[2]==1 else '-'})"
    return dict(ref=ref,gene="intergenic",product=desc,context="intergenic",codon=None,changes="")
