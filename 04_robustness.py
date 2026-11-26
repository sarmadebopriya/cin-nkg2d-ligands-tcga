"""Robustness checks. Everything here is reported, including failures.

  1. Other instability measures (fraction genome altered, CIN70) in place of aneuploidy
  2. Top quarter instead of top third for the Figure 2B split
  3. Purity-high tumors only (top half by purity within each type)
  4. Leave-one-cancer-type-out for the pooled summary
  5. Excluding MSI-high tumors (MANTIS score above config.MSI_HIGH_MANTIS), if an MSI column exists

Run:  python 04_robustness.py
"""
import numpy as np
import pandas as pd

import analysis as A
import config as C


def summarise_fig1(t, analysis, variant, ref=None):
    m = t[~t["is_subset"]]
    sig = m["q_adj"] < C.FDR_ALPHA
    row = dict(analysis=analysis, variant=variant, n_types=len(m),
               median_estimate=m["r_adj"].median(), n_positive=int((m["r_adj"] > 0).sum()),
               n_fdr_sig=int(sig.sum()), n_fdr_sig_positive=int((m.loc[sig, "r_adj"] > 0).sum()))
    if ref is not None:
        j = m.set_index("group")["r_adj"].to_frame("v").join(ref.set_index("group")["r_adj"].to_frame("ref"), how="inner")
        row["same_sign_as_primary"] = int((np.sign(j["v"]) == np.sign(j["ref"])).sum())
        row["n_compared"] = len(j)
    return row


def summarise_int(t, analysis, variant):
    m = t[~t["is_subset"]]
    sig = m["q_int"] < C.FDR_ALPHA
    return dict(analysis=analysis, variant=variant, n_types=len(m),
                median_estimate=m["beta_int"].median(), n_positive=int((m["beta_int"] > 0).sum()),
                n_fdr_sig=int(sig.sum()), n_fdr_sig_positive=int((m.loc[sig, "beta_int"] > 0).sum()))


def main():
    df = A.load_table()
    C.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    # primary reference
    ref = A.figure1_table(df, C.PRIMARY_INSTABILITY)
    if ref.empty:
        raise SystemExit("No cancer type reached MIN_N complete cases.")
    rows.append(summarise_fig1(ref, "fig1 partial rho", f"primary: {C.PRIMARY_INSTABILITY}"))
    rows.append(summarise_int(A.interaction_table(df, C.PRIMARY_INSTABILITY),
                              "fig2A interaction", f"primary: {C.PRIMARY_INSTABILITY}"))

    # 1. other instability measures
    for inst in C.SECONDARY_INSTABILITY:
        if inst not in df.columns or df[inst].notna().sum() == 0:
            rows.append(dict(analysis="fig1 partial rho", variant=f"{inst}: column unavailable"))
            continue
        t = A.figure1_table(df, inst)
        if not t.empty:
            rows.append(summarise_fig1(t, "fig1 partial rho", inst, ref))
        ti = A.interaction_table(df, inst)
        if not ti.empty:
            rows.append(summarise_int(ti, "fig2A interaction", inst))

    # 2. top quarter for the Figure 2B split
    genes = [g for g in C.SENSING_EVASION + C.NK_MARKERS if f"z_{g}" in df.columns]
    for frac, name in [(C.HIGH_FRACTION, "top third (primary)"), (0.25, "top quarter")]:
        se = A.split_effects(df, C.PRIMARY_INSTABILITY, genes, high_fraction=frac)
        if se.empty:
            continue
        m = se[~se["is_subset"]]
        st = m[m["gene"] == "STING1"]
        rows.append(dict(analysis="fig2B split, STING1 effect", variant=name,
                         n_types=len(st), median_estimate=st["effect"].median(),
                         n_positive=int((st["effect"] > 0).sum()),
                         n_fdr_sig=int((m["q"] < C.FDR_ALPHA).sum()),
                         n_fdr_sig_positive=int(((m["q"] < C.FDR_ALPHA) & (m["effect"] > 0)).sum())))

    # 5. excluding MSI-high tumors (near-diploid and immune-rich; can mimic a negative link)
    if "msi" in df.columns and df["msi"].notna().any():
        keep = ~(df["msi"] > C.MSI_HIGH_MANTIS)               # unknown MSI status is kept
        n_drop = int((~keep).sum())
        dfm = df[keep]
        label = f"excluding MSI-high (MANTIS > {C.MSI_HIGH_MANTIS}; {n_drop} tumors)"
        tm = A.figure1_table(dfm, C.PRIMARY_INSTABILITY)
        if not tm.empty:
            rows.append(summarise_fig1(tm, "fig1 partial rho", label, ref))
            j = (ref[["group", "r_adj", "q_adj"]].merge(tm[["group", "n", "r_adj", "q_adj"]], on="group",
                 how="outer", suffixes=("_primary", "_noMSI")))
            j.to_csv(C.RESULTS_DIR / "robustness_msi_comparison.csv", index=False)
            print("Effect of removing MSI-high tumors (purity-adjusted rho):")
            print(j[j["group"].isin(["COAD", "STAD", "UCEC", "READ"])].round(3).to_string(index=False))
        tim = A.interaction_table(dfm, C.PRIMARY_INSTABILITY)
        if not tim.empty:
            rows.append(summarise_int(tim, "fig2A interaction", label))
    else:
        rows.append(dict(analysis="fig1 partial rho", variant="MSI exclusion: column unavailable"))

    # 3. purity-high tumors only
    cut = df.groupby("cancer_type")["purity"].transform(lambda x: x.quantile(C.PURITY_HIGH_QUANTILE))
    dfp = df[df["purity"] >= cut]
    half_n = max(40, C.MIN_N // 2)
    tp = A.figure1_table(dfp, C.PRIMARY_INSTABILITY, min_n=half_n)
    if not tp.empty:
        rows.append(summarise_fig1(tp, "fig1 partial rho", f"purity-high only (n >= {half_n})", ref))
    tpi = A.interaction_table(dfp, C.PRIMARY_INSTABILITY, min_n=half_n)
    if not tpi.empty:
        rows.append(summarise_int(tpi, "fig2A interaction", f"purity-high only (n >= {half_n})"))

    out = pd.DataFrame(rows)
    out.to_csv(C.RESULTS_DIR / "robustness.csv", index=False)

    # 4. leave-one-type-out on the pooled median
    m = ref[~ref["is_subset"]]
    loo = []
    for g in m["group"]:
        rest = m[m["group"] != g]
        loo.append(dict(left_out=g, median_r_adj=rest["r_adj"].median(),
                        n_positive=int((rest["r_adj"] > 0).sum()), n_types=len(rest)))
    loo = pd.DataFrame(loo)
    loo.to_csv(C.RESULTS_DIR / "robustness_leave_one_out.csv", index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print(out.round(3).to_string(index=False))
    print(f"\nLeave-one-type-out median partial rho: min {loo['median_r_adj'].min():.3f}, "
          f"max {loo['median_r_adj'].max():.3f} (all types: {m['r_adj'].median():.3f})")
    flip = loo[np.sign(loo["median_r_adj"]) != np.sign(m["r_adj"].median())]
    if len(flip):
        print(f"WARNING: the sign of the pooled median flips when you drop: {flip['left_out'].tolist()}")
    print("\nWrite down anything that fails in the README. A failure is a result.")


if __name__ == "__main__":
    main()
