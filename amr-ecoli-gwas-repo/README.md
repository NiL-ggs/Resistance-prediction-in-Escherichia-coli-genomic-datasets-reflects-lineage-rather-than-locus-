# Population-structure-aware GWAS of antibiotic resistance in clinical *Escherichia coli*

Code and results for a secondary genome-wide association analysis of the Giessen *E. coli*
dataset (Ren et al., 2022, *Bioinformatics* 38:325–334), covering resistance to ciprofloxacin,
cefotaxime, ceftazidime and gentamicin in 809 isolates at 16,435 SNP positions.

Everything needed to reproduce every number, table and figure in the manuscript and its
supplementary material is in this repository. The two external inputs — the SNP matrix and the
reference annotation — are downloaded by `scripts/00_get_data.sh` rather than redistributed here.
The repository holds the analysis only: the Word documents themselves are written and revised
outside it, and this code produces the numbers, tables and figures they report.

---

## What the analysis does

1. **Deduplicates the public SNP matrix.** It contains 60,936 columns but only 16,446 unique
   positions; all 15,005 multi-column groups are exact row-for-row duplicates. Without this step
   the multiple-testing burden is inflated about 3.7-fold.
2. **Compares predictive accuracy under random and lineage-blocked cross-validation** for four
   classifiers, to separate transferable genotype–phenotype signal from recognition of clonal
   background.
3. **Tests association at four levels of correction for population structure:** none
   (chi-square), linear regression on 2–50 principal components, and a linear mixed model with a
   genome-wide kinship matrix (pyseer). Genomic inflation (λ) is reported for each.
4. **Generates a null distribution** by permuting phenotypes within phylogenetically defined
   lineages (200 permutations per antibiotic) and repeating the mixed-model analysis.
5. **Reconstructs a maximum-likelihood phylogeny** and uses it for Fitch-parsimony homoplasy
   counts and within-lineage (Cochran–Mantel–Haenszel) testing.
6. **Annotates every leading position** against *E. coli* K-12 MG1655 (NC_000913.3), including
   codon number and synonymous/non-synonymous status.
7. **Tests robustness**: inserts a synthetic locus of known effect size into the real genotype
   matrix to confirm the pipeline detects locus effects when present, and repeats the permutation
   null and blocked cross-validation at three lineage-clustering thresholds, with a diagnostic for
   whether each permutation scheme preserves the phenotype's clonal distribution.
8. **Runs positive controls** at known resistance-determining positions, across 27 chromosomal
   genes and regulatory regions with established roles in resistance to the four antibiotics,
   and examines the mixed-model-significant positions outside them as candidate loci.
9. **Separates SNP-level information from population structure** by comparing prediction from
   SNPs alone, from principal components alone and from both, with the components computed
   inside each training fold.

Headline results: classification accuracy falls by 0.11–0.30 AUC for every classifier and every
antibiotic when whole lineages are held out; the number of significant positions falls from
1,788–13,688 (naive) to 7–18 per antibiotic under the mixed model, but those counts are no larger
than a within-lineage permutation null produces (empirical p = 0.73–1.00); and the *gyrA* QRDR,
*parC* and the *ampC* promoter are absent from the public matrix because of its 20 %
variant-frequency filter, so no positive control is available.

---

## Requirements

| Software | Version used | Notes |
|---|---|---|
| Python | 3.12 | |
| pandas | 3.0.2 | |
| numpy | 2.4.4 | |
| scipy | 1.17.1 | |
| scikit-learn | 1.8.0 | PCA |
| statsmodels | 0.15.0 | logistic refit |
| biopython | 1.88 | annotation |
| scikit-bio | 0.7.4 | tree parsing, patristic distances |
| matplotlib | 3.10.8 | |
| openpyxl | 3.1.5 | supplementary workbook |
| pyseer | 1.4.3 | mixed-model association |
| FastTree | 2.2 | phylogeny (`FastTree` or `FastTreeMP` on `PATH`) |

```bash
pip install -r requirements.txt
```

pyseer pulls in optional dependencies (`glmnet_python`, `pybedtools`) that are needed only for
its elastic-net mode, which this analysis does not use. If they fail to build, install pyseer
with `--no-deps` and add `tqdm` separately.

---

## Running it

```bash
bash run_all.sh          # from the repository root, not from scripts/
```

Roughly 6 hours end to end on one core, dominated by the 800 permutation runs (~4.5 h) and the phylogeny (~25 min). Steps 15 and 20 add about 30 minutes. Steps are
independent and numbered, so any one can be rerun on its own provided earlier outputs exist.

| Step | Script | Produces |
|---|---|---|
| 00 | `00_get_data.sh` | `giessen/`, `reference/` |
| 01 | `01_dedup_and_pc_gwas.py` | `results_v2/*_association_results_v2.csv`, `summary_v2.csv` |
| 02 | `02_prepare_matrix.py` | `data.npz`, `snps.fa`, `pyseer/{var.Rtab,pheno.tsv,kinship.tsv}` |
| 03 | `03_pc_sensitivity.py` | `pc_sensitivity.csv`, `PCs50.npy` |
| 04 | `04_pyseer_lmm.sh` | `pyseer/lmm_*.tsv`, `pyseer/pat_*.txt` |
| 05 | `05_build_tree.sh` | `tree.nwk` |
| 06 | `06_homoplasy.py` | `homoplasy.csv` |
| 07 | `07_define_lineages.py` | `tipdist.npy` |
| 08 | `08_within_lineage_cmh.py` | `cmh_lineage.csv`, `lineage_0.1.npy` |
| 09 | `09_logistic_top_hits.py` | `logistic_top50.csv` |
| 10 | `10_make_tables.py` | `out/GWAS_supplementary_tables.xlsx`, `tables.json` |
| 11 | `11_make_figures.py` | `out/Fig_QQ.png`, `Fig_PCA.png`, `Fig_PC_sensitivity.png`, `Fig_coresistance.png` |
| 12 | `12_make_tree_figure.py` | `out/Fig_tree_rings.png` |
| 15 | `15_prediction.py` | `prediction_auc.csv` |
| 16 | `16_permutation_setup.py` | `perm/pheno_perm.tsv` (200 permutations per drug) |
| 17 | `17_permutation_run.sh` | `perm/out/*.tsv` (~4.5 h single-core; re-run to regenerate any incomplete output) |
| 18 | `18_permutation_summary.py` | `perm_summary.csv`, `perm_strongest.csv` |
| 20 | `20_classifier_comparison.py <DRUG>` | `clf_<DRUG>.csv`; concatenate to `classifier_comparison.csv` |
| 21 | `21_make_null_figures.py` | `out/Fig_classifiers.png`, `out/Fig_permutation_null.png` |
| 22 | `22_spike_control.py` + pyseer on `spike/` | `spike/lmm_OR*.tsv` (synthetic-locus positive control) |
| 23 | `23_threshold_sensitivity.py` | `threshold_cv.csv` (blocked CV at 3 thresholds) |
| 24 | `24_threshold_summary.py` | `threshold_perm.csv` |
| 25 | `25_clonality_diagnostic.py` | `clonality_diagnostic.csv` |
| 26 | `26_robustness_figures.py` | `out/Fig_spike.png`, `out/Fig_threshold.png` |
| 02b | `02b_pyseer_inputs.py` | `pyseer/var.Rtab`, `pyseer/pheno.tsv`, `pyseer/kinship.tsv` (also written by step 02; kept separate so steps 04 and 29b can be rerun alone) |
| 27 | `27_snp_vs_pc_prediction.py [n_pcs] [classifiers]` | `snp_vs_pc_auc_*.csv`, `snp_vs_pc_delta_*.csv` → `revision/feature_set_auc.csv`, `revision/feature_set_delta_ci.csv` |
| 28 | `28_known_resistance_genes.py` | `known_gene_positions.csv`, `known_gene_summary.csv`, `known_gene_enrichment.csv` |
| 29a | `29a_gene_level_permutations.py` | `perm_gene/pheno_perm.tsv` (20 within-lineage permutations per drug) |
| 29b | `29b_run_perm.sh` | `perm_gene/lmm_<DRUG>_p<k>.tsv` (~3 h single-core) |
| 29c | `29c_candidate_loci.py` | `candidate_loci.csv`, `candidate_genes.csv`, `gene_level_null.csv` |
| 30 | `30_workflow_figure.py` | `Fig1_workflow.png` |
| 31 | `31_encoding_sensitivity.py` | `encoding_sensitivity.csv` |
| 32 | `32_update_tables_workbook.py` | adds the revision sheets to `GWAS_supplementary_tables.xlsx` |

Steps 27 and 31 are the long ones. Both accept `DRUGS=CIP[,CTX,...]` in the environment, and
step 27 also `MODELS=` (a comma-separated subset of `SNP`, `PC`, `SNP+PC`, `PC (global)`,
`SNP+PC (global)`) and `TAG=` (a suffix for the output file name), so a run can be split per
antibiotic or per feature set and the parts concatenated. For example:

```bash
for d in CIP CTX CTZ GEN; do
  DRUGS=$d MODELS="PC,SNP+PC" TAG=_$d python3 scripts/27_snp_vs_pc_prediction.py 20 RF
done
```
Two points worth knowing before rerunning step 27. Principal components are computed inside each
training fold and the held-out isolates projected onto them; the `PC (global)` and
`SNP+PC (global)` models instead use components computed once from all 809 isolates, which is
optimistic under lineage-blocked evaluation and is reported only for comparison with the earlier
version of this analysis. And in the combined model the standardized component block is rescaled
to the same total variance as the SNP block, so that 20 columns among 16,435 are not swamped by
the shared L2 penalty; without that rescaling the combined model is indistinguishable from the
SNP model by construction.

Gene coordinates for steps 28 and 29c come from the GenBank record when `scripts/00_get_data.sh`
has fetched it, and otherwise from the fallback table in `reference/`, which with the
redistributed `reference/NC_000913.3.fasta` lets those steps run without network access. The
fallback reproduces every gene, codon and amino-acid assignment in Supplementary Tables S4 and
S5; see `reference/README_coordinates.md`.

`scripts/annotate_lib.py` is a helper loaded by steps 09 and 10, and `scripts/annot_mg1655.py`
one loaded by steps 28 and 29c; neither is run directly.

`results/figures/Fig3_manhattan.png` is a main-text figure carried over from the earlier version
of the analysis and is not regenerated by this pipeline.

---

## Inputs

Neither input is redistributed here; both are fetched by step 00.

| File | Source | MD5 |
|---|---|---|
| `cip_ctx_ctz_gen_multi_data.csv` | https://github.com/YunxiaoRen/ML-iAMR | `d768fb095cab002eaa16bde5c3825a77` |
| `cip_ctx_ctz_gen_pheno.csv` | same | `c7a3f28d48c14a5ed62bf444ecbba4dc` |
| `NC_000913.3.gb` | NCBI Nucleotide | `a5fda22bad429251c315ab21e663fa25` |

Check these before interpreting any rerun: the upstream files could change, and the GenBank
record is re-issued periodically, so a different checksum for the reference is not necessarily an
error but should be noted.

### How to read the genotype encoding

The matrix is label-encoded `N=0, A=1, G=2, C=3, T=4`. **N means "no variant called"**, which is
either the MG1655 reference allele or a region absent from or too divergent in that isolate to be
called; the two cannot be separated without the raw reads. For association testing the most
frequent code at each position is the major allele and every other code, N included, is pooled
into one test allele. At the 3,565 multi-allelic positions this pools two distinct states, so the
full allele distribution is reported for all leading positions in Table S3.

---

## Pre-computed results

`results/` holds the outputs of a complete run, so the tables and figures can be checked without
rerunning anything.

- `results/tables/` — per-drug association results (all 16,435 positions, naive and adjusted),
  homoplasy counts, within-lineage p-values, PC-sensitivity, logistic refits, classifier
  comparison, permutation-null summaries, the phylogeny (`tree.nwk`, Newick), and the
  supplementary workbook
- `results/tables/revision/` — results of steps 27–31 (feature-set comparison, encoding
  sensitivity, known resistance genes, candidate loci)
- `results/figures/` — every figure at 300 dpi

---

## Known limitations of what can be verified here

- **Recombination was not tested.** Gubbins and equivalent tools need a whole-genome alignment;
  only the filtered SNP matrix is public. Several reported loci are synonymous changes clustered
  in a narrow interval, a pattern recombination can produce.
- **Plasmid-borne determinants are invisible.** Variants were called against a single chromosomal
  reference, so ESBL and aminoglycoside-modifying genes cannot appear in these results.
- **No positive control exists.** The canonical determinants of all four antibiotics are absent
  from the public matrix; see `results/tables/revision/known_gene_summary.csv` and Supplementary
  Tables S5 and S11. The synthetic-locus spike-in (step 22) is the only demonstration that the
  pipeline detects a locus effect when one is present.
- **No external validation.** An independent collection shares positions with this matrix only
  if its raw reads are called against the same reference with the same 20 % frequency filter.
  Lineage-blocked cross-validation is the closest internal substitute.
- **Permutation validity depends on the lineage threshold.** At a coarse threshold (0.2) the
  permuted phenotype retains only about 60 % of its observed clonal structure, the null collapses
  to zero significant loci, and the comparison becomes uninformative. `25_clonality_diagnostic.py`
  quantifies this; see Section 3.9 of the manuscript. Results at 0.05 and 0.1 agree.
- **Permutation caveat.** Within-lineage permutation preserves relatedness but removes the
  alignment between phenotype and kinship on which the mixed model conditions, making it somewhat
  anti-conservative under this null; the null counts are an upper bound on what structure alone
  produces.

---

## Citation

If you use this code, please cite both the original dataset and this re-analysis:

> Ren Y., Chakraborty T., Doijad S., et al. (2022). Prediction of antimicrobial resistance based
> on whole-genome sequencing and machine learning. *Bioinformatics* 38(2), 325–334.

> [This work — citation to be added on publication.]

## License

Code is released under the MIT License (see `LICENSE`). The Giessen dataset remains under the
terms set by its original authors.