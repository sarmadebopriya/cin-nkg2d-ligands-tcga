# Data (not stored in git)

Put these files in `data/raw/` and `data/`. The full beginner walkthrough, with exact links, is in the project guide. Record the exact dataset names, versions and
download dates in the main README when you fetch them.

| File | What | Where from |
|---|---|---|
| `data/raw/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` | TCGA pan-cancer expression, rows `SYMBOL\|EntrezID` | GDC PanCanAtlas page. A Xena download saved as `data/raw/expression.tsv.gz` also works. |
| `data/raw/PANCAN_ArmCallsAndAneuploidyScore_092817.txt` | Cancer type and aneuploidy score | GDC PanCanAtlas aneuploidy page |
| `data/raw/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | ABSOLUTE tumor purity | GDC PanCanAtlas page |
| `data/raw/mc3.v0.2.8.PUBLIC.maf.gz` (optional) | Mutations, used for TP53 status | GDC PanCanAtlas page |
| `data/raw/brca_clinical.tsv` (optional) | PAM50 subtype for BRCA | cBioPortal, BRCA PanCancer Atlas study clinical data |
| `data/raw/phenotype.tsv` | Built by `python 00_build_phenotype.py` from the files above | made by the script |
| `data/cin70_genes.txt` | CIN70 gene symbols, one per line | Carter et al. 2006, supplementary table. Copy it from the paper, not from memory. |

Only edit `COLS` in `config.py` if you supply your own phenotype file instead of using the builder script.

Cancer type labels must be the same short codes everywhere (BRCA, LUAD, ...). For the BRCA basal-like subset,
`cancer_type` must be `BRCA` and the PAM50 column must contain one of the values in `BRCA_BASAL_VALUES`.

If downloads stall, set `ONLY_TYPES = ["BRCA"]` in `config.py` to run the whole pipeline on breast cancer alone. It still supports the TNBC logic.
