# Step 2: functional annotation and GO term distribution

Integrates the eggNOG-mapper and InterProScan annotations of the Orthopoxvirus
proteomes, resolves redundant GO terms, and tests each conservation group. The
conservation groups come from step 1. The reported result is null, which is why
this step is not titled "enrichment".

eggNOG-mapper v2.1.12 was run in DIAMOND mode, restricted to Viruses, with
E-value <= 1e-5, identity >= 40% and coverage >= 20%, keeping only
non-electronic GO evidence. The test is Fisher's exact test with
Benjamini-Hochberg correction at FDR 0.05. The background population is the set
of 483 proteins carrying at least one GO term, pooled across the four
proteomes.

## The two levels of analysis

The same test runs with two different units of observation.

**Protein level.** One observation per protein. This is exploratory only:
orthologues of the same family are not independent observations, so the test
counts a single ancestral gene once per species and the p-values are
anticonservative. It is reported because it is informative about which terms
the annotation carries, not as evidence.

**Orthogroup level.** One observation per orthogroup. This is the level
reported in the manuscript, and it runs under two aggregation rules:

- "union", where the orthogroup carries the union of the GO terms of its
  annotated members;
- "intersection", where it carries only the terms shared by all of them.

Unannotated members never contribute, so "intersection" is the intersection
across the annotated members. 125 orthogroups have at least one annotated
member, 96 of them in the conserved core and 29 outside it; 87 of the core ones
have exactly four annotated proteins. The union rule tests 850 terms and the
intersection rule 744.

**No term is significant at FDR 0.05 at the orthogroup level**, in any of the
24 comparisons. The protein level does report significant terms; its tables are
labelled as the exploratory level.

The study set is intersected with the population before each test, because an
identifier that carries no annotation would otherwise inflate the study size
with an observation the test cannot evaluate. The run log reports how many
identifiers this drops for every comparison.

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

`../step1/sources/orthogroups.tsv` is read directly, not through `sources/`: it
is what maps each protein to its orthogroup for the orthogroup-level test, with
the two outgroup columns dropped. The R step reads it and
`orthogroups_unassigned_genes.tsv` as well, which is where the
protein-to-species assignment behind the heatmap comes from.

## Requirements

Python 3.11 or newer, with `pandas`, `goatools` and `statsmodels`, at the
versions pinned in the root `requirements.txt`.

The R step needs `dplyr`, `tidyr`, `plotly` and `htmlwidgets`, and `Rscript` on
PATH or `OPXV_RSCRIPT` pointing at the executable.

## Usage

```bash
python run_analysis.py              # the analysis and then the R heatmap
python run_analysis.py pipeline     # the analysis alone
python run_analysis.py enrichment   # the GO test alone, on the existing inputs
```

The test takes no path arguments: its inputs are always
`outputs/goatools/goatools_population.txt`, `outputs/goatools/goatools_study_*.txt`,
`sources/go-basic.obo` and `../step1/sources/orthogroups.tsv`.

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
  the 796 proteins, 483 carry at least one GO term: 333 with a biological
  process, 230 with a cellular component and 323 with a molecular function.
  The table has 796 rows rather than 832 because neither annotation source
  returned anything for 36 proteins of 37 to 86 residues.
- `tables/enhanced_summary_table.tsv`, the per-group summary with GO
  identifiers resolved to term names
- `reports/final_report.txt` and `reports/enhanced_final_report.txt`
- `goatools/goatools_population.txt`, the 483 annotated proteins with their GO
  terms, and one `goatools_study_*.txt` per conservation group, one line per
  protein. The eleven study files together hold the 832 proteins of the four
  proteomes
- `enrichment_dual_level/`, the test results: `protein_level__<group>.tsv`,
  `orthogroup_union__<group>.tsv` and `orthogroup_intersection__<group>.tsv`
  for each conservation group and for the accessory-versus-all comparison,
  plus `summary.tsv` with one row per comparison and `run_log.txt`. Also
  `core_annotated_orthogroups.txt` (96 identifiers) and
  `accessory_annotated_orthogroups.txt` (29), the two partitions the orthogroup
  level is built on
- `visualizations/orthopoxvirus_go_heatmap.html`, the interactive heatmap of
  the top 10 GO terms per category, plus the matrix and top 10 term tables that
  produce it, all written by the R step. The page is standalone, so no
  companion folder travels with it.

Every tested term is written, with its fold enrichment and both p-values, not
only the terms that pass the FDR cutoff. A result table with no significant row
is the reported outcome at the orthogroup level, and it has to be visible as
such rather than inferred from a missing file.
