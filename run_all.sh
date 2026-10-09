#!/usr/bin/env bash
# Full pipeline. Run from the repository root, not from scripts/.
# The analysis only: the manuscript and supplementary Word files are maintained outside
# this repository, which produces the numbers, tables and figures they report.
#
# Runtime is dominated by the phylogeny (step 05, ~25 min), the two permutation runs
# (steps 17 and 29b, several hours single-core) and the feature-set comparison (step 27).
set -euo pipefail
mkdir -p out pyseer perm perm_gene spike

# --- core pipeline -----------------------------------------------------------
bash    scripts/00_get_data.sh              # download SNP matrix and reference
python3 scripts/01_dedup_and_pc_gwas.py     # dedup, naive chi-square, 10-PC linear model
python3 scripts/02_prepare_matrix.py        # data.npz, snps.fa, pyseer inputs
python3 scripts/02b_pyseer_inputs.py        # var.Rtab, pheno.tsv, kinship.tsv (standalone)
python3 scripts/03_pc_sensitivity.py        # lambda and hit counts for 0-50 PCs
bash    scripts/04_pyseer_lmm.sh            # mixed-model association
bash    scripts/05_build_tree.sh            # FastTree phylogeny  (~25 min)
python3 scripts/06_homoplasy.py             # Fitch parsimony, consistency index
python3 scripts/07_define_lineages.py       # patristic distances, lineage clusters
python3 scripts/08_within_lineage_cmh.py    # Cochran-Mantel-Haenszel by lineage
python3 scripts/09_logistic_top_hits.py     # logistic refit, odds ratios
python3 scripts/10_make_tables.py           # supplementary tables workbook
python3 scripts/11_make_figures.py          # QQ, PCA, PC-sensitivity, co-resistance
python3 scripts/12_make_tree_figure.py      # annotated phylogeny
python3 scripts/15_prediction.py            # AUC: LMM loci vs lineage covariates
python3 scripts/16_permutation_setup.py     # within-lineage permuted phenotypes
bash    scripts/17_permutation_run.sh       # 200 permutations x 4 drugs (hours)
python3 scripts/18_permutation_summary.py   # empirical p-values
for d in CIP CTX CTZ GEN; do
  python3 scripts/20_classifier_comparison.py "$d"
done                                         # concatenate clf_*.csv -> classifier_comparison.csv
python3 scripts/21_make_null_figures.py     # classifier and permutation figures
python3 scripts/22_spike_control.py         # synthetic-locus positive control (then pyseer on spike/)
python3 scripts/23_threshold_sensitivity.py # blocked CV at three lineage thresholds
python3 scripts/24_threshold_summary.py     # permutation null at three thresholds
python3 scripts/25_clonality_diagnostic.py  # is the permutation scheme valid at each threshold?
python3 scripts/26_robustness_figures.py    # spike and threshold figures
python3 scripts/27_snp_vs_pc_prediction.py 20 LR   # SNPs vs PCs vs both (see README for DRUGS/MODELS/TAG)
python3 scripts/28_known_resistance_genes.py       # 27 known resistance genes and regions
python3 scripts/29a_gene_level_permutations.py     # 20 within-lineage permutations per drug
bash    scripts/29b_run_perm.sh                    # mixed model on each permuted phenotype (hours)
python3 scripts/29c_candidate_loci.py              # candidate loci outside the known genes
python3 scripts/30_workflow_figure.py              # workflow figure
python3 scripts/31_encoding_sensitivity.py         # binary vs one-hot vs integer codes
python3 scripts/32_update_tables_workbook.py       # add the revision sheets to the workbook

echo "Done. Tables and figures are in out/; pre-computed copies are in results/."
