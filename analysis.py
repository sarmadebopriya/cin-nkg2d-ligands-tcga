"""Shared statistics and helpers. No file paths or choices live here; see config.py."""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

import config as C

# ---------------------------------------------------------------- transforms
def rank_int(x):
    """Rank-based inverse normal transform. NaN stays NaN."""
    s = pd.Series(np.asarray(x, dtype=float))
    ok = s.notna()
    out = pd.Series(np.nan, index=s.index)
    r = s[ok].rank(method="average")
    out[ok] = stats.norm.ppf((r - 0.5) / ok.sum())
    return out.to_numpy()


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


def fdr(p, alpha=C.FDR_ALPHA):
    """Benjamini-Hochberg. NaN p-values stay NaN and are left out of the family."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() > 0:
        q[ok] = multipletests(p[ok], alpha=alpha, method="fdr_bh")[1]
    return q


# ---------------------------------------------------------------- correlation
def _fisher_ci(r, se):
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def spearman_ci(x, y):
    """Spearman r, p, and 95% CI (Bonett-Wright standard error)."""
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    n = len(d)
    if n < 10:
        return dict(r=np.nan, p=np.nan, lo=np.nan, hi=np.nan, n=n)
    r, p = stats.spearmanr(d["x"], d["y"])
    se = np.sqrt((1 + r**2 / 2) / (n - 3))
    lo, hi = _fisher_ci(r, se)
    return dict(r=float(r), p=float(p), lo=lo, hi=hi, n=n)


def partial_spearman(x, y, covars):
    """Spearman partial correlation: rank all variables, remove covariates by
    linear regression on ranks, correlate residuals. CI by Fisher z with n-3-k df."""
    cov = pd.DataFrame(covars)
    d = pd.concat([pd.Series(np.asarray(x), name="x"), pd.Series(np.asarray(y), name="y"),
                   cov.reset_index(drop=True)], axis=1).dropna()
    n, k = len(d), cov.shape[1]
    if n < 10 + k:
        return dict(r=np.nan, p=np.nan, lo=np.nan, hi=np.nan, n=n)
    ranks = d.rank()
    Z = np.column_stack([np.ones(n), ranks[cov.columns].to_numpy()])
    def resid(v):
        b, *_ = np.linalg.lstsq(Z, v, rcond=None)
        return v - Z @ b
    rx, ry = resid(ranks["x"].to_numpy()), resid(ranks["y"].to_numpy())
    r = float(np.corrcoef(rx, ry)[0, 1])
    df = n - 2 - k
    t = r * np.sqrt(df / max(1e-12, 1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    lo, hi = _fisher_ci(r, 1 / np.sqrt(n - 3 - k))
    return dict(r=r, p=p, lo=lo, hi=hi, n=n)


# ---------------------------------------------------------------- group handling
def iter_groups(df, include_subsets=True):
    """Yield (label, sub-dataframe, is_subset). Main groups are cancer types.
    The BRCA basal-like subset is added as a labeled extra, outside the FDR family."""
    for ct, sub in df.groupby("cancer_type", sort=True):
        yield str(ct), sub, False
    if include_subsets and "pam50" in df.columns:
        brca = df[df["cancer_type"] == C.BRCA_TYPE_LABEL]
        basal = brca[brca["pam50"].isin(C.BRCA_BASAL_VALUES)]
        if len(basal) > 0:
            yield f"{C.BRCA_TYPE_LABEL}-basal (subset)", basal, True


def instability_column(name):
    assert name in ("aneuploidy", "fga", "cin70"), name
    return name


# ---------------------------------------------------------------- Figure 1 statistics
def figure1_table(df, instability, outcome="ligand_score", purity_col=None, min_n=C.MIN_N):
    """Per-group Spearman and purity-adjusted partial Spearman, with FDR over main groups only."""
    purity_col = purity_col or "purity"
    icol = instability_column(instability)
    rows, skipped = [], []
    for label, sub, is_subset in iter_groups(df):
        need = [icol, outcome, purity_col]
        s = sub.dropna(subset=need)
        floor = C.MIN_N_SUBSET if is_subset else min_n
        if len(s) < floor:
            skipped.append(f"{label} (n={len(s)})")
            continue
        raw = spearman_ci(s[icol], s[outcome])
        adj = partial_spearman(s[icol], s[outcome], {"purity": s[purity_col]})
        rows.append(dict(group=label, is_subset=is_subset, n=len(s),
                         r=raw["r"], p=raw["p"], lo=raw["lo"], hi=raw["hi"],
                         r_adj=adj["r"], p_adj=adj["p"], lo_adj=adj["lo"], hi_adj=adj["hi"]))
    t = pd.DataFrame(rows)
    t.attrs["skipped"] = skipped
    if t.empty:
        return t
    main = ~t["is_subset"]
    t["q"] = np.nan
    t["q_adj"] = np.nan
    t.loc[main, "q"] = fdr(t.loc[main, "p"])
    t.loc[main, "q_adj"] = fdr(t.loc[main, "p_adj"])
    return t.sort_values("r_adj").reset_index(drop=True)


# ---------------------------------------------------------------- Figure 2 statistics
def interaction_table(df, instability, outcome="ligand_score", min_n=C.MIN_N, sting_col="z_STING1",
                      purity_col=None):
    """Per-group OLS: outcome ~ instability * STING1 + purity (+ TP53). All continuous
    variables are rank-inverse-normal transformed within the group. HC3 robust errors.
    The interaction coefficient asks whether the instability-ligand slope depends on STING1."""
    purity_col = purity_col or "purity"
    icol = instability_column(instability)
    tp53 = "tp53_mut" if "tp53_mut" in df.columns else None
    rows, skipped = [], []
    for label, sub, is_subset in iter_groups(df):
        need = [icol, outcome, purity_col, sting_col] + ([tp53] if tp53 else [])
        s = sub.dropna(subset=need).copy()
        floor = C.MIN_N_SUBSET if is_subset else min_n
        if len(s) < floor:
            skipped.append(f"{label} (n={len(s)})")
            continue
        m = pd.DataFrame({
            "y": rank_int(s[outcome]), "x": rank_int(s[icol]),
            "s": rank_int(s[sting_col]), "purity": zscore(s[purity_col]),
        })
        terms = ["x * s", "purity"]
        used_tp53 = False
        if tp53 and s[tp53].nunique() > 1:
            m["tp53"] = s[tp53].to_numpy(dtype=float)
            terms.append("tp53")
            used_tp53 = True
        fit = smf.ols("y ~ " + " + ".join(terms), data=m).fit(cov_type="HC3")
        rows.append(dict(group=label, is_subset=is_subset, n=len(m), tp53_adjusted=used_tp53,
                         beta_inst=fit.params["x"], beta_sting=fit.params["s"],
                         beta_int=fit.params["x:s"],
                         se_int=fit.bse["x:s"], p_int=fit.pvalues["x:s"],
                         lo_int=fit.conf_int().loc["x:s", 0], hi_int=fit.conf_int().loc["x:s", 1]))
    t = pd.DataFrame(rows)
    t.attrs["skipped"] = skipped
    if t.empty:
        return t
    t["q_int"] = np.nan
    main = ~t["is_subset"]
    t.loc[main, "q_int"] = fdr(t.loc[main, "p_int"])
    return t.sort_values("beta_int").reset_index(drop=True)


def _purity_residual(values, purity):
    """Rank both, regress ranks on purity ranks, return residuals (same idea as partial Spearman)."""
    r = pd.Series(values).rank().to_numpy(dtype=float)
    p = pd.Series(purity).rank().to_numpy(dtype=float)
    Z = np.column_stack([np.ones(len(p)), p])
    b, *_ = np.linalg.lstsq(Z, r, rcond=None)
    return r - Z @ b


def split_effects(df, instability, genes, outcome="ligand_score", high_fraction=C.HIGH_FRACTION,
                  min_group=C.MIN_GROUP, min_n=C.MIN_N, adjust_purity=True):
    """Within each group: take instability-high tumors (top `high_fraction`), split at the
    median ligand score into ligand-high and ligand-low, and compare each gene (z-scored
    expression) with Mann-Whitney. Effect size is the rank-biserial correlation, positive
    when the gene is higher in the ligand-HIGH group. FDR over all main-group x gene tests.
    With adjust_purity=True (default) the ligand score and each gene are first made
    purity-free within the cancer type (rank residuals), so tumor purity cannot create the split."""
    icol = instability_column(instability)
    rows = []
    for label, sub, is_subset in iter_groups(df):
        need = [icol, outcome] + (["purity"] if adjust_purity else [])
        s = sub.dropna(subset=need).copy()
        floor_n = C.MIN_N_SUBSET if is_subset else min_n
        floor_g = C.MIN_GROUP_SUBSET if is_subset else min_group
        if len(s) < floor_n:
            continue
        out_col = outcome
        if adjust_purity:
            s["_out_adj"] = _purity_residual(s[outcome], s["purity"])
            out_col = "_out_adj"
            for g in genes:
                if f"z_{g}" in s.columns:
                    ok = s[f"z_{g}"].notna()
                    s.loc[ok, f"_adj_{g}"] = _purity_residual(s.loc[ok, f"z_{g}"], s.loc[ok, "purity"])
        cut = s[icol].quantile(1 - high_fraction)
        hi = s[s[icol] >= cut]
        med = hi[out_col].median()
        lig_hi, lig_lo = hi[hi[out_col] > med], hi[hi[out_col] <= med]
        if len(lig_hi) < floor_g or len(lig_lo) < floor_g:
            continue
        for g in genes:
            col = f"_adj_{g}" if adjust_purity else f"z_{g}"
            if col not in s.columns:
                continue
            a, b = lig_hi[col].dropna(), lig_lo[col].dropna()
            if len(a) < floor_g or len(b) < floor_g:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(dict(group=label, is_subset=is_subset, gene=g, n_ligand_high=len(a),
                             n_ligand_low=len(b), effect=2 * u / (len(a) * len(b)) - 1, p=p))
    t = pd.DataFrame(rows)
    if t.empty:
        return t
    t["q"] = np.nan
    main = ~t["is_subset"]
    t.loc[main, "q"] = fdr(t.loc[main, "p"])
    return t


def load_table():
    path = C.PROC_DIR / "analysis_table.csv.gz"
    if not path.exists():
        raise SystemExit(f"{path} not found. Run 01_prepare.py first.")
    return pd.read_csv(path)
