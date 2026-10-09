import pandas as pd
import numpy as np
from scipy import stats
import re, json, os
from collections import defaultdict

DATA_DIR = "giessen/Giessen_dataset"
OUT_DIR = "results_v2"
os.makedirs(OUT_DIR, exist_ok=True)

print("Loading data...")
geno = pd.read_csv(f"{DATA_DIR}/cip_ctx_ctz_gen_multi_data.csv").set_index("prename")
pheno = pd.read_csv(f"{DATA_DIR}/cip_ctx_ctz_gen_pheno.csv").set_index("prename")
common = geno.index.intersection(pheno.index)
geno = geno.loc[common]
pheno = pheno.loc[common]
n_samples = len(common)
print(f"Samples: {n_samples}")

# ---------- STEP 1: DEDUPLICATE IDENTICAL COLUMNS ----------
print("\n--- Deduplicating SNP columns ---")
pos_to_cols = defaultdict(list)
for c in geno.columns:
    pos = int(re.match(r"X(\d+)", c).group(1))
    pos_to_cols[pos].append(c)

print(f"Raw columns: {geno.shape[1]}, unique positions: {len(pos_to_cols)}")

# Verify every multi-column group is an exact duplicate before collapsing (safety check)
kept_cols = []
positions = []
n_verified_exact = 0
n_found_different = 0
for pos, cs in pos_to_cols.items():
    if len(cs) > 1:
        sub = geno[cs].values
        if not np.all(sub == sub[:, [0]]):
            n_found_different += 1
            # Keep all distinct columns if they actually differ (shouldn't happen, but don't silently drop data)
            seen = {}
            for c in cs:
                key = tuple(geno[c].values)
                if key not in seen:
                    seen[key] = c
                    kept_cols.append(c)
                    positions.append(pos)
            continue
        else:
            n_verified_exact += 1
    kept_cols.append(cs[0])
    positions.append(pos)

print(f"Confirmed exact duplicates collapsed: {n_verified_exact} positions")
print(f"Positions where duplicates were NOT identical (kept separately): {n_found_different}")
print(f"Final deduplicated column count: {len(kept_cols)}")

G_df = geno[kept_cols].copy()
G_df.columns = [f"pos_{p}" for p in positions]
positions_arr = np.array(positions)
n_snps = G_df.shape[1]
assert n_snps == len(pos_to_cols), "Mismatch after dedup — should equal unique position count"

# ---------- Recode major/minor allele (binary) ----------
G = G_df.values.astype(np.int8)
counts_per_col = np.zeros((5, n_snps), dtype=np.int64)
for code in range(5):
    counts_per_col[code] = (G == code).sum(axis=0)
major_allele = counts_per_col.argmax(axis=0).astype(np.int8)
Gbin = (G != major_allele[np.newaxis, :]).astype(np.float64)

variant_counts = Gbin.sum(axis=0)
polymorphic_mask = (variant_counts > 0) & (variant_counts < n_samples)
Gbin_poly = Gbin[:, polymorphic_mask]
positions_poly = positions_arr[polymorphic_mask]
n_poly = Gbin_poly.shape[1]
print(f"Polymorphic, testable, DEDUPLICATED SNPs: {n_poly}")

drugs = ["CIP", "CTX", "CTZ", "GEN"]
drug_names = {"CIP": "Ciprofloxacin", "CTX": "Cefotaxime", "CTZ": "Ceftazidime", "GEN": "Gentamicin"}

def vectorized_chi2(Gbin_poly, y):
    n = len(y)
    y = y.astype(np.float64)
    y1 = y.sum(); y0 = n - y1
    n1 = Gbin_poly.sum(axis=0); n0 = n - n1
    n11 = y @ Gbin_poly
    n10 = n1 - n11
    n01 = y1 - n11
    n00 = n0 - n01
    e11 = n1*y1/n; e10 = n1*y0/n; e01 = n0*y1/n; e00 = n0*y0/n
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2 = ((n11-e11)**2/e11 + (n10-e10)**2/e10 + (n01-e01)**2/e01 + (n00-e00)**2/e00)
    chi2 = np.nan_to_num(chi2, nan=0.0, posinf=0.0, neginf=0.0)
    pvals = stats.chi2.sf(chi2, df=1)
    a,b,c,d = n11,n10,n01,n00
    zero_cell = (a==0)|(b==0)|(c==0)|(d==0)
    a_,b_,c_,d_ = a.copy(),b.copy(),c.copy(),d.copy()
    a_[zero_cell]+=0.5; b_[zero_cell]+=0.5; c_[zero_cell]+=0.5; d_[zero_cell]+=0.5
    with np.errstate(divide='ignore', invalid='ignore'):
        odds_ratio = (a_*d_)/(b_*c_)
    return chi2, pvals, odds_ratio, n11, n10, n01, n00

def bh_fdr(pvals):
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n)
    out[order] = q
    return out

# ---------- STEP 2: PCA for population-structure correction ----------
print("\n--- Computing PCA on deduplicated genotype matrix (lineage structure) ---")
from sklearn.decomposition import PCA
Gc = Gbin_poly - Gbin_poly.mean(axis=0, keepdims=True)
pca = PCA(n_components=10, random_state=0)
PCs = pca.fit_transform(Gc)
explained = pca.explained_variance_ratio_
print("Variance explained by top 10 PCs:", np.round(explained, 4))
print("Cumulative:", np.round(np.cumsum(explained), 4))

np.save(f"{OUT_DIR}/PCs.npy", PCs)
pd.DataFrame({"PC": range(1, 11), "variance_explained": explained}).to_csv(
    f"{OUT_DIR}/pca_variance_explained.csv", index=False)

# Design matrix with intercept + PCs, for projecting out population structure (Frisch-Waugh-Lovell)
X = np.column_stack([np.ones(n_samples), PCs])
XtX_inv = np.linalg.pinv(X.T @ X)
Hat = X @ XtX_inv @ X.T   # n x n projection (hat) matrix onto covariate space
k_covariates = X.shape[1]  # intercept + 10 PCs = 11

def pc_adjusted_test(Gbin_poly, y, Hat, k_covariates):
    n = Gbin_poly.shape[0]
    y = y.astype(np.float64)
    y_resid = y - Hat @ y
    G_resid = Gbin_poly - Hat @ Gbin_poly   # residualize every SNP column in one matmul

    Gss = np.sum(G_resid * G_resid, axis=0)
    valid = Gss > 1e-8
    beta = np.zeros(Gbin_poly.shape[1])
    beta[valid] = (G_resid[:, valid].T @ y_resid) / Gss[valid]

    resid_full = y_resid[:, None] - beta[None, :] * G_resid  # residuals of full model, per SNP
    sse = np.sum(resid_full**2, axis=0)
    df_resid = n - k_covariates - 1
    var_beta = np.full(Gbin_poly.shape[1], np.inf)
    var_beta[valid] = sse[valid] / df_resid / Gss[valid]
    se_beta = np.sqrt(var_beta)

    with np.errstate(divide='ignore', invalid='ignore'):
        tstat = beta / se_beta
    tstat = np.nan_to_num(tstat, nan=0.0)
    pvals = 2 * stats.t.sf(np.abs(tstat), df=df_resid)
    pvals[~valid] = 1.0
    return beta, tstat, pvals

summary_rows = []
for drug in drugs:
    y = pheno.loc[common, drug].values.astype(np.int8)
    n_resistant = int(y.sum())
    print(f"\n=== {drug} ({drug_names[drug]}) ===")

    # --- Naive (uncorrected) test, on CLEAN deduplicated data ---
    chi2, pvals, odds_ratio, n11, n10, n01, n00 = vectorized_chi2(Gbin_poly, y)
    qvals = bh_fdr(pvals)
    bonferroni_thresh = 0.05 / n_poly
    n_bonf = int((pvals < bonferroni_thresh).sum())
    n_fdr = int((qvals < 0.05).sum())
    print(f"  Naive (dedup):        Bonferroni-sig={n_bonf}, FDR5%-sig={n_fdr}  (out of {n_poly} unique SNPs)")

    # --- PC-adjusted (population-structure-corrected) test ---
    beta, tstat, pvals_adj = pc_adjusted_test(Gbin_poly, y, Hat, k_covariates)
    qvals_adj = bh_fdr(pvals_adj)
    n_bonf_adj = int((pvals_adj < bonferroni_thresh).sum())
    n_fdr_adj = int((qvals_adj < 0.05).sum())
    print(f"  PC-adjusted:          Bonferroni-sig={n_bonf_adj}, FDR5%-sig={n_fdr_adj}  (out of {n_poly} unique SNPs)")

    df_out = pd.DataFrame({
        "position": positions_poly,
        "chi2_naive": chi2, "p_naive": pvals, "q_naive": qvals, "odds_ratio_naive": odds_ratio,
        "beta_pc_adjusted": beta, "t_pc_adjusted": tstat,
        "p_pc_adjusted": pvals_adj, "q_pc_adjusted": qvals_adj,
    }).sort_values("p_pc_adjusted").reset_index(drop=True)
    df_out.to_csv(f"{OUT_DIR}/{drug}_association_results_v2.csv", index=False)

    summary_rows.append({
        "drug": drug, "drug_name": drug_names[drug],
        "n_resistant": n_resistant, "n_susceptible": n_samples - n_resistant,
        "n_unique_snps_tested": n_poly,
        "naive_bonferroni_sig": n_bonf, "naive_fdr5pct_sig": n_fdr,
        "pc_adjusted_bonferroni_sig": n_bonf_adj, "pc_adjusted_fdr5pct_sig": n_fdr_adj,
        "reduction_factor_fdr": round(n_fdr / n_fdr_adj, 2) if n_fdr_adj > 0 else float('inf'),
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{OUT_DIR}/summary_v2.csv", index=False)
print("\n=== FINAL SUMMARY (deduplicated + PC-adjusted) ===")
print(summary_df.to_string(index=False))
