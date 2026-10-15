# Step 2: functional annotation and GO enrichment

Integrates the eggNOG-mapper and InterProScan annotations of the Orthopoxvirus
proteomes, resolves redundant GO terms, and runs enrichment per conservation
group. The conservation groups come from step 1.

eggNOG-mapper v2.1.12 was run in DIAMOND mode, restricted to Viruses, with
E-value <= 1e-5, identity >= 40% and coverage >= 20%, keeping only
non-electronic GO evidence. Enrichment uses Fisher's exact test with
Benjamini-Hochberg correction at FDR 0.05, with the proteome of each virus as
the background population.

## Inputs

Required in `sources/`:

- `eggnog_annotations.tsv`
- `interpro_orthopoxvirus.tsv`
- `classified_protein_ids.txt` (the step 1 report, under its short name)
- `go-basic.obo`
- `all_refseq_proteins.fasta`
- `list_*.txt`, one file per conservation group. The group and the species are
  read from the file name, so these cannot be renamed on their own.

Optional:

- `ReactomePathways.txt`, used to resolve pathway identifiers to names

The R step also reads `orthogroups.tsv` and `orthogroups_unassigned_genes.tsv`
from `../step1/sources/`, which is where the protein-to-species assignment
behind the heatmap comes from.

## Requirements

Python 3.11 or newer, with `pandas` and `goatools`.

The R step needs `Rscript` with `dplyr`, `tidyr`, `plotly` and `htmlwidgets`.

## Usage

```bash
python run_analysis.py            
```

## Outputs

Under `outputs/`:

- `tables/step3_final_consolidated_annotations.tsv`, the consolidated
  annotation table, and the intermediate tables that produce it. One row per
  protein, with each source in its own column: `eggNOG_Description`,
  `Interpro_Accession` and `Interpro_Description`, then the GO terms split by
  namespace into `GO_BP`, `GO_CC` and `GO_MF`. Each GO cell carries the
  identifier and the term name together, as `GO:0005576 (extracellular
  region)`, and the terms it lists are already pruned of parents that a more
  specific descendant makes redundant. A namespace with no term reads `-`. Of
  the 796 proteins, 241 carry at least one GO term: 150 with a biological
  process, 168 with a cellular component and 59 with a molecular function.
- `tables/enhanced_summary_table.tsv`, the per-group summary with GO
  identifiers resolved to term names
- `reports/final_report.txt` and `reports/enhanced_final_report.txt`
- `goatools/goatools_population.txt` and one `goatools_study_*.txt` per
  conservation group
- `visualizations/orthopoxvirus_go_heatmap.html`, the interactive heatmap of
  the top 10 GO terms per category, plus the matrix and top 10 term tables that
  produce it, all written by the R step. The page is standalone, so no
  companion folder travels with it.

Enrichment writes a result file only when at least one term passes the FDR
cutoff, so a missing file for a group means nothing was significant for it.
