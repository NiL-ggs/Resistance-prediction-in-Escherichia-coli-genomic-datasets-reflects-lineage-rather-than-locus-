"""Helper: annotate matrix positions against NC_000913.3 using reference/mg1655_genes.json
(named CDS), reference/mg1655_genemark_cds.json (all predicted CDS) and the reference
sequence. Used by steps 28-29."""
import json, numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
CODE={0:"N",1:"A",2:"G",3:"C",4:"T"}
import os
def _from_genbank(path):
    """Named CDS coordinates from the GenBank record, when it has been downloaded (step 00)."""
    rec=next(SeqIO.parse(path,"genbank")); out={}
    for f in rec.features:
        if f.type=="CDS" and "gene" in f.qualifiers:
            out.setdefault(f.qualifiers["gene"][0],dict(start=int(f.location.start)+1,end=int(f.location.end),
                           strand="+" if f.location.strand==1 else "-"))
    return out
GB="reference/NC_000913.3.gb"; FA="reference/NC_000913.3.fasta"
_REC=None
if os.path.exists(GB) and os.path.getsize(GB)>1_000_000:
    _REC=next(SeqIO.parse(GB,"genbank"))
    GENES=_from_genbank(GB); SOURCE="GenBank NC_000913.3"
else:
    # fallback: coordinate table reconstructed and validated as described in
    # reference/README_coordinates.md, used when the GenBank record is unavailable
    GENES=json.load(open("reference/mg1655_genes.json")); SOURCE="reconstructed table"
    if not os.path.exists(FA):
        raise SystemExit(
            "No reference sequence found. Run scripts/00_get_data.sh to fetch "
            "reference/NC_000913.3.gb (preferred), or place the NC_000913.3 sequence at "
            + FA + " — see reference/README_coordinates.md.")
GENES.update({k:v for k,v in {"gyrB":dict(start=3877705,end=3880119,strand="-"),"acrR":dict(start=485761,end=486408,strand="+"),
              "acrA":dict(start=484426,end=485619,strand="-"),"marA":dict(start=1619574,end=1619957,strand="+"),
              "soxS":dict(start=4277060,end=4277383,strand="-"),"soxR":dict(start=4277469,end=4277933,strand="+"),
              "fusA":dict(start=3471400,end=3473514,strand="-")}.items() if k not in GENES})
REF=str(_REC.seq) if _REC is not None else str(next(SeqIO.parse(FA,"fasta")).seq)
_names=sorted(GENES.items(),key=lambda kv:kv[1]["start"])
_starts=np.array([v["start"] for _,v in _names])
def gene_at(p):
    hits=[k for k,v in _names if v["start"]<=p<=v["end"]]
    return hits[0] if hits else None
_GM=json.load(open("reference/mg1655_genemark_cds.json")) if SOURCE!="GenBank NC_000913.3" else {}
def flanks(p):
    i=np.searchsorted(_starts,p)
    left=_names[i-1][0] if i>0 else None; right=_names[i][0] if i<len(_names) else None
    if any(s<=p<=e for s,e,_ in _GM.values()): return f"unnamed CDS between {left} and {right}"
    return f"intergenic {left}–{right}"
def effect(p,gname,states):
    g=GENES[gname]
    off=p-g["start"] if g["strand"]=="+" else g["end"]-p
    cn=off//3+1; k=off%3
    cs=g["start"]+3*(cn-1) if g["strand"]=="+" else g["end"]-3*cn+1
    cod=REF[cs-1:cs+2]; cod=cod if g["strand"]=="+" else str(Seq(cod).reverse_complement())
    a0=str(Seq(cod).translate()); out=[]
    for code in sorted(set(states)-{0}):
        b=CODE[code] if g["strand"]=="+" else str(Seq(CODE[code]).complement())
        alt=cod[:k]+b+cod[k+1:]
        if alt==cod: continue
        a1=str(Seq(alt).translate())
        out.append(f"{CODE[code]}:{a0}{cn}{a1}"+(" (syn)" if a1==a0 else " (non-syn)"))
    return "; ".join(out) if out else "reference state only"
def annotate(p,states):
    g=gene_at(p)
    if g: return g,effect(p,g,states)
    return flanks(p),"intergenic or unnamed CDS"
