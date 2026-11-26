"""Which cancer type is significant in the instability x STING1 interaction after MSI-high tumors are removed?
Saves results/robustness_msi_fig2a.csv and prints the types with q below 0.05."""
import config as C
import analysis as A

df = A.load_table()
if "msi" not in df.columns or not df["msi"].notna().any():
    raise SystemExit("No MSI column in the analysis table.")
dfm = df[~(df["msi"] > C.MSI_HIGH_MANTIS)]
prim = A.interaction_table(df, C.PRIMARY_INSTABILITY)
nom = A.interaction_table(dfm, C.PRIMARY_INSTABILITY)
j = prim[["group", "beta_int", "q_int"]].merge(
    nom[["group", "n", "beta_int", "lo_int", "hi_int", "q_int"]],
    on="group", how="outer", suffixes=("_primary", "_noMSI")).sort_values("q_int_noMSI")
C.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
j.to_csv(C.RESULTS_DIR / "robustness_msi_fig2a.csv", index=False)
print("Interaction after removing MSI-high tumors, sorted by q:")
print(j.head(8).round(3).to_string(index=False))
print("\nSignificant (q below %.2f):" % C.FDR_ALPHA)
print(j[j["q_int_noMSI"] < C.FDR_ALPHA].round(3).to_string(index=False))
