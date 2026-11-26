"""Load TCGA expression + phenotype, filter to one primary tumor per patient, z-score within
each cancer type, compute CIN70 and ligand scores, and write one analysis table.

Run:  python 01_prepare.py
"""
import sys
import numpy as np
import pandas as pd

import config as C

report = []
def say(msg=""):
    print(msg)
    report.append(str(msg))


def resolve_genes():
    """Map each wanted gene to the symbol present in the expression file."""
    wanted = list(dict.fromkeys(C.LIGANDS + C.SENSING_EVASION + C.NK_MARKERS))
    candidates = {g: C.ALIASES.get(g, [g]) for g in wanted}
    return wanted, candidates


def read_expression(genes_needed):
    """Read the expression matrix in chunks, keeping only the rows we need."""
    if not C.EXPR_FILE.exists():
        sys.exit(f"Expression file not found: {C.EXPR_FILE}\nSee data/README.md for how to fetch it.")
    print(f"Reading expression from {C.EXPR_FILE.name} (this can take a few minutes)...")
    keep = []
    for chunk in pd.read_csv(C.EXPR_FILE, sep="\t", index_col=0, chunksize=2000):
        # GDC pan-cancer files label rows "SYMBOL|EntrezID"; keep only the symbol
        parts = chunk.index.astype(str).str.split("|")
        sym = parts.str[0].str.upper()
        ent = parts.str[1] if parts.str.len().max() > 1 else pd.Series("", index=chunk.index)
        chunk.index = [C.ENTREZ_TO_GENE.get(str(e), s) for s, e in zip(sym, ent)]
        hit = chunk[chunk.index.isin(genes_needed)]
        if len(hit):
            keep.append(hit)
    if not keep:
        sys.exit("None of the requested genes were found. Is the file genes-by-samples with gene symbols in column 1?")
    expr = pd.concat(keep)
    expr = expr[~expr.index.duplicated(keep="first")]
    return expr.T  # samples x genes


def sample_type(barcode):
    return barcode[13:15] if len(barcode) >= 15 else None


def main():
    C.PROC_DIR.mkdir(parents=True, exist_ok=True)
    C.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ---- CIN70 list (must come from the paper, not memory)
    if C.CIN70_FILE.exists():
        cin70 = [g.strip().upper() for g in C.CIN70_FILE.read_text().split() if g.strip()]
        say(f"CIN70 list: {len(cin70)} genes read from {C.CIN70_FILE.name}")
    else:
        cin70 = []
        say(f"No CIN70 list at {C.CIN70_FILE.name}: skipping the CIN70 robustness check. "
            "(Optional: add the Carter 2006 gene list, one symbol per line, and rerun.)")

    wanted, candidates = resolve_genes()
    all_names = {n.upper() for names in candidates.values() for n in names} | set(cin70)
    expr = read_expression(all_names)
    present = set(expr.columns)

    # ---- resolve gene symbols
    resolved, missing = {}, []
    for g in wanted:
        hit = next((n.upper() for n in candidates[g] if n.upper() in present), None)
        if hit:
            resolved[g] = hit
        else:
            missing.append(g)
    say("Gene check:")
    for g in wanted:
        say(f"  {g:8s} -> {resolved.get(g, 'MISSING')}")
    if [g for g in missing if g in C.REQUIRED_GENES]:
        sys.exit(f"Required genes missing: {[g for g in missing if g in C.REQUIRED_GENES]}\n"
                 "Send this printout to Claude. Do not edit the script.")
    if missing:
        say(f"WARNING: not found, will be skipped: {missing}")
    cin_present = [g for g in cin70 if g in present]
    if cin70:
        say(f"CIN70 genes found in expression file: {len(cin_present)} of {len(cin70)}")
        if len(cin_present) < 50:
            sys.exit("Fewer than 50 CIN70 genes found. Check symbols or the gene list.")

    # ---- log scale check on the genes we use
    use_cols = list(resolved.values()) + cin_present
    vals = expr[use_cols].to_numpy(dtype=float)
    if np.nanmax(vals) > 50:
        say("Expression looks linear (max > 50): applying log2(x + 1).")
        expr[use_cols] = np.log2(expr[use_cols].astype(float) + 1)
    else:
        say("Expression looks log-scale already (max <= 50): no transform applied.")

    # ---- phenotype and sample ids
    if not C.PHENO_FILE.exists():
        sys.exit(f"Phenotype file not found: {C.PHENO_FILE}")
    ph = pd.read_csv(C.PHENO_FILE, sep="\t")
    rename = {v: k for k, v in C.COLS.items() if v is not None and v in ph.columns}
    absent = [k for k, v in C.COLS.items() if v is not None and v not in ph.columns]
    ph = ph.rename(columns=rename)
    for need in ("id", "cancer_type", "purity", "aneuploidy"):
        if need not in ph.columns:
            sys.exit(f"Phenotype file lacks the '{need}' column (config.COLS['{need}']={C.COLS[need]!r}).")
    optional_absent = [k for k in absent if k in ("fga", "tp53_mut", "pam50", "msi")]
    if optional_absent:
        say(f"WARNING: optional phenotype columns not found, related analyses will be skipped: {optional_absent}")

    expr = expr.copy()
    expr.index.name = "sample_full"
    expr = expr.reset_index()
    expr["sample"] = expr["sample_full"].str[:15]
    expr["patient"] = expr["sample_full"].str[:12]
    expr["stype"] = expr["sample_full"].map(sample_type)

    if C.PHENO_KEY == "sample":
        ph["sample"] = ph["id"].astype(str).str[:15]
        df = expr.merge(ph.drop(columns=["id"]), on="sample", how="inner")
    else:
        ph["patient"] = ph["id"].astype(str).str[:12]
        df = expr.merge(ph.drop(columns=["id"]), on="patient", how="inner")
    say(f"Samples in expression: {len(expr)}; matched to phenotype: {len(df)}")
    if len(df) < 0.5 * len(expr):
        say("WARNING: fewer than half of the expression samples matched the phenotype file. "
            "Check that both use TCGA barcodes (PHENO_KEY in config.py).")
    if C.ONLY_TYPES:
        df = df[df["cancer_type"].isin(C.ONLY_TYPES)]
        say(f"ONLY_TYPES is set: keeping {C.ONLY_TYPES}; {len(df)} samples left")

    # ---- keep primary tumors only, one per patient
    def keep_priority(row):
        order = C.SAMPLE_TYPE_PRIORITY.get(row["cancer_type"], C.SAMPLE_TYPE_PRIORITY["default"])
        return order.index(row["stype"]) if row["stype"] in order else np.nan
    df["prio"] = df.apply(keep_priority, axis=1)
    n_before = len(df)
    df = df.dropna(subset=["prio"]).sort_values(["patient", "prio"]).drop_duplicates("patient", keep="first")
    say(f"After keeping primary tumors (SKCM metastatic) and one sample per patient: {len(df)} of {n_before}")

    # ---- z-score within cancer type, then scores
    gene_cols = {g: s for g, s in resolved.items()}
    for g, sym in gene_cols.items():
        df[f"z_{g}"] = df.groupby("cancer_type")[sym].transform(
            lambda x: (x - x.mean()) / x.std(ddof=1))
    if cin_present:
        zc = df.groupby("cancer_type")[cin_present].transform(lambda x: (x - x.mean()) / x.std(ddof=1))
        df["cin70"] = zc.mean(axis=1)
    lig_z = [f"z_{g}" for g in C.LIGANDS if g in gene_cols]
    df["ligand_score"] = df[lig_z].mean(axis=1)
    say(f"Ligand score = mean z of: {[g for g in C.LIGANDS if g in gene_cols]}")

    # ---- output table and counts
    keep = ["sample", "patient", "cancer_type", "purity", "aneuploidy", "fga", "tp53_mut", "pam50", "msi",
            "cin70", "ligand_score"] + [f"z_{g}" for g in gene_cols]
    keep = [c for c in keep if c in df.columns]
    out = df[keep].copy()
    out.to_csv(C.PROC_DIR / "analysis_table.csv.gz", index=False)

    counts = out.groupby("cancer_type").agg(
        n=("sample", "size"),
        n_purity=("purity", lambda x: x.notna().sum()),
        n_aneuploidy=("aneuploidy", lambda x: x.notna().sum()),
    )
    counts["n_complete"] = out.dropna(subset=["purity", "aneuploidy"]).groupby("cancer_type").size()
    counts = counts.fillna(0).astype(int)
    counts["analysed"] = counts["n_complete"] >= C.MIN_N
    counts.to_csv(C.RESULTS_DIR / "sample_counts.csv")
    say("\nTumors per cancer type (analysed = at least MIN_N complete cases):")
    say(counts.to_string())
    say(f"\nTypes analysed: {int(counts['analysed'].sum())} of {len(counts)} (MIN_N = {C.MIN_N})")

    (C.RESULTS_DIR / "prepare_report.txt").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
