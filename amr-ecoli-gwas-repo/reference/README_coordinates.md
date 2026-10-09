# MG1655 (NC_000913.3) gene coordinates used by the revision analyses

`scripts/annot_mg1655.py` reads the GenBank record (`reference/NC_000913.3.gb`, downloaded by
`scripts/00_get_data.sh`) when it is present, and the files in this directory otherwise. The
fallback exists only so that the gene-set analyses added in revision can be rerun in an
environment without access to NCBI.

- `mg1655_genes.json` — named CDS coordinates (start, end, strand). Reconstructed from
  locus-tagged coding-variant annotation on NC_000913.3: for each variant, the genomic
  position and the coding-sequence offset in the HGVS `c.` notation fix the gene's start (plus
  strand) or end (minus strand), and the annotated CDS length fixes the other end. Genes with a
  single annotated variant were resolved by requiring a start and a stop codon in the reference
  sequence in exactly one of the two orientations.
- `mg1655_genemark_cds.json` — all predicted CDS on the reference (GeneMark.hmm, E. coli K-12
  model), used to flag matrix positions that fall in an unnamed coding sequence and to supply
  seven genes absent from the source above (gyrB, acrA, acrR, marA, soxR, soxS, fusA). Each of
  those seven was confirmed by translating the predicted CDS and matching its N-terminal
  protein sequence.
- `NC_000913.3.fasta` — the NC_000913.3 sequence (public domain, from NCBI), used for codon
  extraction when the GenBank record has not been downloaded. It is the one input redistributed
  here, so that the fallback path works without network access.

Validation: these coordinates reproduce every gene, codon and amino-acid assignment in
Supplementary Tables S4 and S5, which were produced from the GenBank record (for example gyrA
codons 85-783, yhgF codon 521, mgtA codon 610, ampC D367A, marR G103S and Y137H).
