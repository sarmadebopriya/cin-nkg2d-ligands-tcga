"""Figure 2: among unstable tumors, what separates ligand-high from ligand-low?

Panel A (primary): per cancer type, does STING1 expression change how instability relates to
ligand score? OLS with an instability x STING1 interaction, adjusted for purity (and TP53 if
available). Variables are rank-inverse-normal transformed; HC3 robust errors; FDR over types.

Panel B (descriptive): within instability-high tumors, split at the median ligand score and
compare CGAS, STING1, TREX1, ADAM10, ADAM17 and NK markers. The ligand score and every gene are
first made purity-free within each cancer type. NK markers mainly track immune infiltrate, so
read them as descriptive only. An unadjusted version is saved as a CSV for comparison.

Run:  python 03_figure2.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import analysis as A
import config as C
import plotting as P


def main():
    df = A.load_table()
    C.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    inst = C.PRIMARY_INSTABILITY

    it = A.interaction_table(df, inst, C.PRIMARY_OUTCOME)
    if it.empty:
        raise SystemExit("No cancer type reached MIN_N complete cases.")
    it.to_csv(C.RESULTS_DIR / f"fig2a_interaction_{inst}.csv", index=False)

    genes = [g for g in C.SENSING_EVASION + C.NK_MARKERS if f"z_{g}" in df.columns]
    se = A.split_effects(df, inst, genes, C.PRIMARY_OUTCOME)           # purity-adjusted (primary)
    if not se.empty:
        se.to_csv(C.RESULTS_DIR / f"fig2b_split_{inst}.csv", index=False)
    se_raw = A.split_effects(df, inst, genes, C.PRIMARY_OUTCOME, adjust_purity=False)  # for comparison only
    if not se_raw.empty:
        se_raw.to_csv(C.RESULTS_DIR / f"fig2b_split_{inst}_unadjusted.csv", index=False)

    ncol = 2 if not se.empty else 1
    fig = plt.figure(figsize=(6 + 10 * (ncol - 1), 0.28 * len(it) + 2.2))
    ax = fig.add_subplot(1, ncol, 1)
    P.forest(ax, it["group"].tolist(), it["beta_int"].to_numpy(), it["lo_int"].to_numpy(),
             it["hi_int"].to_numpy(), it["q_int"].to_numpy(), it["is_subset"].tolist(),
             "A. Instability x STING1 interaction", "interaction coefficient (95% CI)")
    ax.text(0.5, -0.12 - 0.5 / max(len(it), 1), "positive: ligand rises faster with instability when STING1 is high",
            transform=ax.transAxes, fontsize=7, color=P.GREY, ha="center", va="top")

    if not se.empty:
        ax2 = fig.add_subplot(1, 2, 2)
        piv = se.pivot(index="gene", columns="group", values="effect").reindex(genes)
        qv = se.pivot(index="gene", columns="group", values="q").reindex(genes)
        annot = piv.round(2).astype(str)
        annot = annot.where(~(qv < C.FDR_ALPHA), annot + "*")
        annot = annot.where(piv.notna(), "")
        sns.heatmap(piv, annot=annot, fmt="", cmap="vlag", center=0, vmin=-0.8, vmax=0.8,
                    cbar_kws={"label": "rank-biserial effect\n(+ = higher in ligand-high)"}, ax=ax2,
                    annot_kws={"fontsize": 6})
        ax2.set_title(f"B. Instability-high tumors: ligand-high vs ligand-low (purity-adjusted)\n* FDR < {C.FDR_ALPHA}; "
                      "subset column has no FDR", fontsize=9)
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        plt.setp(ax2.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    P.save(fig, "figure2_sting_status_and_ligand_loss")

    main_rows = it[~it["is_subset"]]
    if it.attrs.get("skipped"):
        print(f"Skipped for low n: {it.attrs['skipped']}")
    print(f"Interaction (instability x STING1), {len(main_rows)} cancer types: "
          f"median beta {main_rows['beta_int'].median():.3f}; "
          f"FDR<{C.FDR_ALPHA}: {(main_rows['q_int'] < C.FDR_ALPHA).sum()} "
          f"({(main_rows[main_rows['q_int'] < C.FDR_ALPHA]['beta_int'] > 0).sum()} positive)")
    if not se.empty:
        s = se[(~se["is_subset"]) & (se["q"] < C.FDR_ALPHA)]
        print(f"Split comparison (purity-adjusted): {len(s)} gene x type tests at FDR<{C.FDR_ALPHA}")
    if not se_raw.empty:
        sr = se_raw[(~se_raw["is_subset"]) & (se_raw["q"] < C.FDR_ALPHA)]
        print(f"Split comparison (NOT purity-adjusted, for comparison): {len(sr)} tests at FDR<{C.FDR_ALPHA}")
    print("Saved results/ and figures/ for Figure 2.")


if __name__ == "__main__":
    main()
