"""Figure 1: does genomic instability track NKG2D ligand expression within each cancer type?

Primary: aneuploidy score vs ligand score, Spearman, per cancer type, with a purity-adjusted
partial Spearman beside it. FDR over cancer types. Per-gene results go to the supplement.

Run:  python 02_figure1.py
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

    t = A.figure1_table(df, inst, C.PRIMARY_OUTCOME)
    if t.empty:
        raise SystemExit("No cancer type reached MIN_N complete cases.")
    t.to_csv(C.RESULTS_DIR / f"fig1_{inst}_ligand_score.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(9, 0.28 * len(t) + 1.8), sharey=True)
    lab = t["group"].tolist()
    sub = t["is_subset"].tolist()
    P.forest(axes[0], lab, t["r"].to_numpy(), t["lo"].to_numpy(), t["hi"].to_numpy(),
             t["q"].to_numpy(), sub, "Spearman", "rho (95% CI)")
    P.forest(axes[1], lab, t["r_adj"].to_numpy(), t["lo_adj"].to_numpy(), t["hi_adj"].to_numpy(),
             t["q_adj"].to_numpy(), sub, "Partial Spearman, purity-adjusted", "rho (95% CI)")
    fig.suptitle(f"Instability ({inst}) vs ligand score, per cancer type", fontsize=11, y=1.0)
    fig.text(0.5, -0.02, f"Filled = FDR < {C.FDR_ALPHA}; open = not significant; "
             "square (orange) = labeled subset, outside FDR", ha="center", fontsize=8, color=P.GREY)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    P.save(fig, "figure1_instability_vs_ligand")

    # ---- supplement: per-gene, purity-adjusted, one cell per cancer type x ligand
    rows = []
    for g in C.LIGANDS:
        col = f"z_{g}"
        if col not in df.columns:
            continue
        tg = A.figure1_table(df, inst, outcome=col)
        if tg.empty:
            continue
        tg["gene"] = g
        rows.append(tg)
    if rows:
        pg = pd.concat(rows)
        pg.to_csv(C.RESULTS_DIR / f"fig1_{inst}_per_gene.csv", index=False)
        piv = pg.pivot(index="group", columns="gene", values="r_adj")
        qv = pg.pivot(index="group", columns="gene", values="q_adj")
        annot = piv.round(2).astype(str)
        annot = annot.where(~(qv < C.FDR_ALPHA), annot + "*")
        fig, ax = plt.subplots(figsize=(0.9 * piv.shape[1] + 2, 0.3 * piv.shape[0] + 1.5))
        sns.heatmap(piv, annot=annot, fmt="", cmap="vlag", center=0, vmin=-0.6, vmax=0.6,
                    cbar_kws={"label": "partial Spearman rho"}, ax=ax, annot_kws={"fontsize": 7})
        ax.set_title(f"Supplement: per-ligand, purity-adjusted ({inst}). * FDR < {C.FDR_ALPHA}", fontsize=10)
        ax.set_ylabel("")
        P.save(fig, "figureS1_per_ligand")

    main_rows = t[~t["is_subset"]]
    print(f"Cancer types analysed: {len(main_rows)}")
    if t.attrs.get("skipped"):
        print(f"Skipped for low n: {t.attrs['skipped']}")
    print(f"Purity-adjusted rho, median {main_rows['r_adj'].median():.3f}; "
          f"positive in {(main_rows['r_adj'] > 0).sum()}; "
          f"FDR<{C.FDR_ALPHA}: {(main_rows['q_adj'] < C.FDR_ALPHA).sum()} "
          f"({(main_rows[main_rows['q_adj'] < C.FDR_ALPHA]['r_adj'] > 0).sum()} positive)")
    print(f"Saved results/ and figures/ for Figure 1.")


if __name__ == "__main__":
    main()
