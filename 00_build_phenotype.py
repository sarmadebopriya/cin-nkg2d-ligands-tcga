"""Build data/raw/phenotype.tsv from the PanCanAtlas files you downloaded.

Needed:  PANCAN_ArmCallsAndAneuploidyScore_092817.txt   (cancer type + aneuploidy score)
         TCGA_mastercalls.abs_tables_JSedit.fixed.txt   (tumor purity)
Optional: mc3.v0.2.8.PUBLIC.maf.gz  (TP53 mutation status; big, slow)
          brca_clinical.tsv         (PAM50 subtype for the BRCA basal-like subset)

Run:  python 00_build_phenotype.py
If a column cannot be found the script prints every column name it saw. Nothing is guessed silently.
"""
import re
import sys

import pandas as pd

import config as C

NON_SILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
              "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation",
              "Translation_Start_Site"}


def norm(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def find_col(df, path, what, exact=(), contains=()):
    """Return the first column whose normalized name equals one of `exact`
    or contains one of `contains`. Otherwise stop and show the columns."""
    for c in df.columns:                      # exact names win over partial matches
        if norm(c) in exact:
            return c
    for c in df.columns:
        if any(k in norm(c) for k in contains):
            return c
    cols = list(df.columns)
    sys.exit(f"\nCould not find the {what} column in {path.name}.\n"
             f"Columns seen ({len(cols)}): {cols[:40]}{' ...' if len(cols) > 40 else ''}\n"
             "Open the file, find the right column name, and tell me (or edit this script).")


def read_table(path, what):
    if not path.exists():
        sys.exit(f"Missing file for {what}: {path}\nSee GUIDE step 'Download the data'.")
    return pd.read_csv(path, sep="\t", low_memory=False)


def main():
    # ---- aneuploidy score and cancer type
    an = read_table(C.ANEUPLOIDY_FILE, "aneuploidy score")
    id_a = an.columns[0]
    type_a = find_col(an, C.ANEUPLOIDY_FILE, "cancer type", exact=("type", "cancertype", "disease"))
    score_a = find_col(an, C.ANEUPLOIDY_FILE, "aneuploidy score", contains=("aneuploidyscore",))
    ph = pd.DataFrame({
        "sample": an[id_a].astype(str).str[:15],
        "cancer_type": an[type_a].astype(str),
        "aneuploidy_score": pd.to_numeric(an[score_a], errors="coerce"),
    }).drop_duplicates("sample")
    print(f"Aneuploidy file: {len(ph)} samples; columns used: {id_a!r}, {type_a!r}, {score_a!r}")

    # ---- purity
    pu = read_table(C.PURITY_FILE, "tumor purity")
    id_p = find_col(pu, C.PURITY_FILE, "sample ID", exact=("sample",), contains=("array",))
    pur_p = find_col(pu, C.PURITY_FILE, "purity", exact=("purity",))
    pu = pd.DataFrame({"sample": pu[id_p].astype(str).str[:15],
                       "purity": pd.to_numeric(pu[pur_p], errors="coerce")})
    pu = pu.groupby("sample", as_index=False)["purity"].mean()
    ph = ph.merge(pu, on="sample", how="left")
    print(f"Purity file: {len(pu)} samples; columns used: {id_p!r}, {pur_p!r}")

    # ---- optional: TP53 from the MC3 MAF
    if C.MAF_FILE.exists():
        print("Reading the MC3 mutation file in chunks (several minutes)...")
        seen, mutated = set(), set()
        for ch in pd.read_csv(C.MAF_FILE, sep="\t", comment="#", low_memory=False, chunksize=500_000,
                              usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"]):
            ids = ch["Tumor_Sample_Barcode"].astype(str).str[:15]
            seen.update(ids.unique())
            hit = ch[(ch["Hugo_Symbol"] == "TP53") & ch["Variant_Classification"].isin(NON_SILENT)]
            mutated.update(hit["Tumor_Sample_Barcode"].astype(str).str[:15].unique())
        ph["tp53_mut"] = [1.0 if s in mutated else (0.0 if s in seen else float("nan")) for s in ph["sample"]]
        print(f"TP53: {int((ph['tp53_mut'] == 1).sum())} mutated, {int((ph['tp53_mut'] == 0).sum())} wild type, "
              f"{int(ph['tp53_mut'].isna().sum())} unknown")
    else:
        print(f"No MAF file at {C.MAF_FILE.name}: skipping TP53 (the model then runs without it).")

    # ---- optional: PAM50 for BRCA
    if C.PAM50_FILE.exists():
        cl = read_table(C.PAM50_FILE, "PAM50 subtype")
        sub_c = find_col(cl, C.PAM50_FILE, "subtype", exact=("subtype", "pam50", "pam50subtype"))
        label = cl[sub_c].astype(str).str.replace("BRCA_", "", regex=False)
        # Prefer a Sample ID column; fall back to Patient ID (cBioPortal may omit Sample ID).
        samp_c = next((c for c in cl.columns if norm(c) in ("sampleid", "sample", "samplebarcode")), None)
        pat_c = next((c for c in cl.columns if norm(c) in ("patientid", "patient", "bcrpatientbarcode")), None)
        if samp_c is not None:
            sub = pd.DataFrame({"sample": cl[samp_c].astype(str).str[:15], "pam50": label})
            sub = sub.drop_duplicates("sample")
            ph = ph.merge(sub, on="sample", how="left")
            print(f"PAM50: matched on {samp_c!r}")
        elif pat_c is not None:
            sub = pd.DataFrame({"_pat": cl[pat_c].astype(str).str[:12], "pam50": label}).drop_duplicates("_pat")
            ph["_pat"] = ph["sample"].str[:12]
            ph = ph.merge(sub, on="_pat", how="left").drop(columns="_pat")
            print(f"PAM50: no Sample ID column, matched on {pat_c!r} (first 12 characters of the barcode)")
        else:
            sys.exit(f"\nNo Sample ID or Patient ID column in {C.PAM50_FILE.name}.\n"
                     f"Columns seen: {list(cl.columns)}\nRe-download with Patient ID or Sample ID ticked.")
        n_basal = int(((ph["pam50"] == "Basal") & (ph["cancer_type"] == C.BRCA_TYPE_LABEL)).sum())
        print(f"PAM50: {n_basal} BRCA basal-like tumors matched")
        if n_basal == 0:
            print("WARNING: no basal-like tumors matched. Check the sample ID format and subtype labels.")
    else:
        print(f"No PAM50 file at {C.PAM50_FILE.name}: the BRCA basal-like subset will be skipped.")

    # ---- optional: fraction genome altered and MSI score from the cBioPortal pan-cancer table
    if C.PANCAN_CLINICAL_FILE.exists():
        pc = read_table(C.PANCAN_CLINICAL_FILE, "pan-cancer clinical table")
        fga_c = find_col(pc, C.PANCAN_CLINICAL_FILE, "fraction genome altered",
                         exact=("fractiongenomealtered", "fga"))
        msi_c = next((c for c in pc.columns if norm(c) in ("msimantisscore", "mantis")), None)
        pcs = pd.DataFrame({"fga": pd.to_numeric(pc[fga_c], errors="coerce")})
        if msi_c is not None:
            pcs["msi_mantis"] = pd.to_numeric(pc[msi_c], errors="coerce")
        samp_c = next((c for c in pc.columns if norm(c) in ("sampleid", "sample", "samplebarcode")), None)
        pat_c = next((c for c in pc.columns if norm(c) in ("patientid", "patient")), None)
        if samp_c is not None:
            pcs["sample"] = pc[samp_c].astype(str).str[:15]
            pcs = pcs.drop_duplicates("sample")
            ph = ph.merge(pcs, on="sample", how="left")
            how_matched = repr(samp_c)
        elif pat_c is not None:
            pcs["_pat"] = pc[pat_c].astype(str).str[:12]
            pcs = pcs.drop_duplicates("_pat")
            ph["_pat"] = ph["sample"].str[:12]
            ph = ph.merge(pcs, on="_pat", how="left").drop(columns="_pat")
            how_matched = f"{pat_c!r} (first 12 characters)"
        else:
            sys.exit(f"\nNo Sample ID or Patient ID column in {C.PANCAN_CLINICAL_FILE.name}.\n"
                     f"Columns seen: {list(pc.columns)}")
        print(f"Pan-cancer clinical table matched on {how_matched}: FGA for {int(ph['fga'].notna().sum())} samples"
              + (f", MSI MANTIS for {int(ph['msi_mantis'].notna().sum())}" if "msi_mantis" in ph else ", no MSI MANTIS column found"))
    else:
        print(f"No pan-cancer clinical table at {C.PANCAN_CLINICAL_FILE.name}: skipping FGA and MSI checks.")

    C.RAW_DIR.mkdir(parents=True, exist_ok=True)
    out = C.RAW_DIR / "phenotype.tsv"
    ph.to_csv(out, sep="\t", index=False)
    print(f"\nWrote {out}  ({len(ph)} rows)")
    print("Rows per cancer type (top):")
    print(ph["cancer_type"].value_counts().head(40).to_string())
    print(f"Missing purity: {int(ph['purity'].isna().sum())}; missing aneuploidy: {int(ph['aneuploidy_score'].isna().sum())}")


if __name__ == "__main__":
    main()
