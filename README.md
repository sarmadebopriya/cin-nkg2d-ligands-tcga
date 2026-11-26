# Genomic instability and NKG2D ligand expression across TCGA tumors

**Status:** pipeline written and tested on synthetic data only. No real result is reported here yet.
Fill the bracketed parts after the real run.

## Question
Within each cancer type, do tumors with more genomic instability express more NKG2D ligand mRNA, and does
STING1 expression change that relationship? This is a descriptive check of the premise behind my PhD aim on
stress-state control of surface NKG2D ligands in triple-negative breast cancer. It cannot test causes.

## Data
GDC PanCanAtlas supplemental files (record the download date for each):
- RNA: EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv (batch-corrected pan-cancer RSEM)
- Aneuploidy scores: PANCAN_ArmCallsAndAneuploidyScore_092817.txt (Taylor et al. 2018)
- Purity: TCGA_mastercalls.abs_tables_JSedit.fixed.txt (ABSOLUTE)
- TP53 status (optional): mc3.v0.2.8.PUBLIC.maf.gz
- PAM50 (optional): [source and date]
- CIN70 list: Carter et al. 2006, supplementary table

## Choices fixed before looking at results (config.py)
- Primary instability measure: aneuploidy score. Secondary: fraction of genome altered (only if you add an
  `fga` column; the default downloads do not include it). CIN70 is a further check. CIN70 is largely a
  proliferation signature, so it is only a check.
- Primary outcome: ligand score, the mean within-type z-score of MICA, MICB, ULBP1, ULBP2, ULBP3,
  RAET1E (ULBP4) and RAET1G (ULBP5). Per-ligand results are supplementary.
- Every analysis runs within a cancer type, never pooled. Types need at least 150 complete tumors.
- "Instability-high" is the top third within each type. Ligand-high vs ligand-low splits at the median
  ligand score inside the instability-high group.
- Purity is a covariate everywhere. TP53 mutation is a covariate in the interaction model.
- Multiple testing: Benjamini-Hochberg over cancer types. The BRCA basal-like subset is labeled and sits
  outside the FDR family.

## Run order
```
pip install -r requirements.txt
python 00_build_phenotype.py   # merge cancer type, aneuploidy score, purity (+ TP53, PAM50 if downloaded)
python 01_prepare.py      # filter, z-score, score, write data_processed/analysis_table.csv.gz
python 02_figure1.py      # instability vs ligand score per cancer type (+ per-ligand supplement)
python 03_figure2.py      # instability x STING1 interaction (A) and ligand-high vs ligand-low (B)
python 04_robustness.py   # other instability measures, top quarter, purity-high, leave-one-type-out
```
To test the code without real data: `python tests/make_synthetic.py`, then run the four scripts with
`CIN_DATA_DIR=tests/synthetic_data CIN_OUT_DIR=tests/synthetic_out` (the synthetic files mimic the real
download formats, so `00_build_phenotype.py` runs on them too). Synthetic output must never be reported.

## Methods in brief
Primary tumors only, one sample per patient (SKCM: metastatic samples preferred). Expression z-scored within
each type. Figure 1: Spearman correlation with a Bonett-Wright confidence interval, and a purity-adjusted
partial Spearman correlation (ranks, residualized on purity). Figure 2A: per type, OLS of ligand score on
instability, STING1, their product, purity and TP53, all continuous variables rank-inverse-normal transformed,
HC3 robust errors; the product term asks whether STING1 changes the instability-ligand slope.
Figure 2B: Mann-Whitney tests with rank-biserial effect sizes.

## Results
[Fill in after the real run. State exactly what you found, including weak or null results.]
- Figure 1: [n types, median partial rho, how many positive, how many pass FDR]
- Figure 2A: [interaction coefficient direction and how many types pass FDR]
- Figure 2B: [descriptive only; NK markers mostly track immune infiltrate]
- BRCA basal-like subset: [report separately, labeled as a subset]
- Robustness: [what held, what did not]

## What this cannot resolve
- Bulk mRNA is not surface protein, and TCGA cannot see shedding. ADAM10 and ADAM17 mRNA say little about
  sheddase activity. Ligand regulation is also largely post-transcriptional, so weak correlations are plausible.
- Aneuploidy and fraction of genome altered measure chromosomal instability. They are not replication stress or
  cytosolic DNA, which is what the proposed aim is about.
- STING1, CGAS and the NK markers come heavily from immune and stromal cells in bulk tumors. Purity
  adjustment helps but does not isolate tumor-intrinsic expression.
- Microsatellite-unstable tumors have low aneuploidy but strong immune infiltrate, which can distort
  patterns in colon, uterine and stomach cancers.
- PAM50 basal-like is a proxy for TNBC, not the same set of tumors.
- All results are associations in bulk tumors. None of them show a cause.
