"""Build a SYNTHETIC dataset that mimics the file formats of the real downloads, to test the pipeline.
The numbers are random with a planted signal. Never report results from this data.

Run:  python tests/make_synthetic.py
Then (Mac/Linux):
  export CIN_DATA_DIR=tests/synthetic_data CIN_OUT_DIR=tests/synthetic_out
  python 00_build_phenotype.py && python 01_prepare.py && python 02_figure1.py && python 03_figure2.py && python 04_robustness.py
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd

rng = np.random.default_rng(1)
out = Path(__file__).resolve().parent / "synthetic_data"
(out / "raw").mkdir(parents=True, exist_ok=True)

types = {"BRCA": 420, "LUAD": 300, "COAD": 260, "KIRC": 220, "LGG": 200, "SKCM": 180,
         "OV": 170, "HNSC": 160, "STAD": 155, "GBM": 90, "CHOL": 40}
ligands = ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G"]
other = ["CGAS", "STING1", "TREX1", "ADAM10", "ADAM17", "NCR1", "KLRD1", "PRF1", "GZMB"]
cin = [f"CIN{i}" for i in range(1, 71)]          # placeholder names for the test only
filler = [f"GENE{i}" for i in range(300)]
# Real GDC files use older symbols for two genes; mimic that to test the alias logic
old_symbol = {"STING1": "TMEM173", "CGAS": "MB21D1"}

cols, aneu, absolute, maf, clin = {}, [], [], [], []
pid = 0
for t, n in types.items():
    base = rng.normal(0, 1, len(ligands) + len(other)) * 0.3 + 5
    for _ in range(n):
        pid += 1
        patient = f"TCGA-AA-{pid:04d}"
        purity = float(np.clip(rng.beta(5, 2), 0.1, 0.99))
        inst = rng.normal(0, 1)
        tp53 = int(rng.random() < 0.3 + 0.2 * (inst > 0.5))
        sting = rng.normal(0, 1) + 0.5 * (1 - purity)
        infil = rng.normal(0, 1) + 1.5 * (1 - purity)
        lig = 0.25 * inst + 0.25 * inst * sting + 0.1 * tp53 + rng.normal(0, 1)
        g = {}
        for k, gname in enumerate(ligands):
            g[gname] = base[k] + 0.8 * lig + rng.normal(0, 1)
        g["CGAS"] = 5 + 0.3 * infil + rng.normal(0, 1)
        g["STING1"] = 5 + sting
        for gname in ["TREX1", "ADAM10", "ADAM17"]:
            g[gname] = 5 + rng.normal(0, 1)
        for gname in ["NCR1", "KLRD1", "PRF1", "GZMB"]:
            g[gname] = 5 + 0.8 * infil + rng.normal(0, 1)
        for c in cin:
            g[c] = 6 + 0.6 * inst + 0.5 * rng.normal(0, 1)
        for f in filler:
            g[f] = 5 + rng.normal(0, 1)
        stype = "06" if (t == "SKCM" and rng.random() < 0.8) else "01"
        sample15 = f"{patient}-{stype}"
        cols[f"{sample15}A-11R-A29S-07"] = g            # full aliquot barcode like the GDC file
        if rng.random() < 0.04:                          # a normal sample and a second tumor vial
            cols[f"{patient}-11A-11R-A29S-07"] = g
            cols[f"{patient}-{stype}B-11R-A29S-07"] = g
        aneu.append({"Sample": sample15, "Type": t,
                     "AneuploidyScore(AS)": int(np.clip(round(10 + 4 * inst + rng.normal(0, 2)), 0, 39))})
        absolute.append({"array": f"{sample15}A", "sample": sample15,
                         "purity": purity if rng.random() > 0.03 else np.nan})
        if rng.random() > 0.02:
            maf.append({"Hugo_Symbol": "TP53" if tp53 else "KRAS",
                        "Variant_Classification": "Missense_Mutation", "Tumor_Sample_Barcode": f"{sample15}A-11D-A1B2-09"})
            maf.append({"Hugo_Symbol": "TP53", "Variant_Classification": "Silent",
                        "Tumor_Sample_Barcode": f"{sample15}A-11D-A1B2-09"}) if rng.random() < 0.1 else None
        if t == "BRCA":
            clin.append({"Sample ID": sample15, "Subtype": rng.choice(
                ["BRCA_Basal", "BRCA_LumA", "BRCA_LumB", "BRCA_Her2", "BRCA_Normal"], p=[.2, .4, .2, .1, .1])})

expr = pd.DataFrame(cols)
expr.index = [f"{old_symbol.get(g, g)}|{1000 + i}" for i, g in enumerate(expr.index)]
expr = 2 ** expr - 1                                   # real file is linear scale, not log
expr.index.name = "gene_id"
expr.to_csv(out / "raw" / "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", sep="\t")
aneuploidy = pd.DataFrame(aneu)
aneuploidy["Type2"] = aneuploidy["Type"]               # real file repeats a Type column
aneuploidy.rename(columns={"Type2": "Type"}).to_csv(
    out / "raw" / "PANCAN_ArmCallsAndAneuploidyScore_092817.txt", sep="\t", index=False)
pd.DataFrame(absolute).to_csv(out / "raw" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t", index=False)
pd.DataFrame([m for m in maf if m]).to_csv(out / "raw" / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", index=False,
                                           compression="gzip")
pd.DataFrame(clin).to_csv(out / "raw" / "brca_clinical.tsv", sep="\t", index=False)
(out / "cin70_genes.txt").write_text("\n".join(cin) + "\n")
print("synthetic data written to", out, "| expression", expr.shape)
