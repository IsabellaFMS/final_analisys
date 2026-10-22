# Step 2: functional annotation and GO term distribution

Integrates the eggNOG-mapper and InterProScan annotations of the Orthopoxvirus
proteomes, resolves redundant GO terms, and tests each conservation group. The
conservation groups come from step 1.

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

The study set is intersected with the population before each test, because an
identifier that carries no annotation would otherwise inflate the study size
with an observation the test cannot evaluate.

## Inputs

Required in `sources/`:

- `eggnog_annotations.tsv`
- `interpro_orthopoxvirus.tsv`
- `go-basic.obo`
- `all_refseq_proteins.fasta`
- `list_*.txt`, one file per conservation group. The group and the species are
  read from the file name, so these cannot be renamed on their own.

Optional:

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

- `tables/final_consolidated_annotations.tsv`,
- `tables/enhanced_summary_table.tsv`,
- `reports/final_report.txt` and `reports/enhanced_final_report.txt`
- `goatools/goatools_population.txt`, the 483 annotated proteins with their GO
  terms, and one `goatools_study_*.txt` per conservation group, one line per
  protein.
- `enrichment_dual_level/`, the test results: `protein_level__<group>.tsv`,
  `orthogroup_union__<group>.tsv` and `orthogroup_intersection__<group>.tsv`
  for each conservation group and for the accessory-versus-all comparison,
  plus `summary.tsv` with one row per comparison and `run_log.txt`.
- `visualizations/orthopoxvirus_go_heatmap.html`. The page is standalone, so no
  companion folder travels with it.
