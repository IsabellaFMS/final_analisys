# Step 1: orthogroup conservation map

Consolidates the OrthoFinder orthogroup tables and classifies every viral
protein by how many of the four orthopoxviruses share its orthogroup. The
classification report produced here is the input that steps 2, 4, 5 and 6 use
to label their own results.

## Inputs

Both files come from the same OrthoFinder run:

- `sources/orthogroups.tsv`
- `sources/orthogroups_unassigned_genes.tsv`

The two files share 128 `Orthogroup` IDs. The loader groups rows by
`Orthogroup`, merges the protein lists per species and keeps one record per
group, so those 128 groups are not counted twice.

## Excluded species

`Humanalphaherpesvirus1_proteins` and `swinepox_proteins` are present in the
OrthoFinder output but are dropped before any metric is computed. They are not
orthopoxviruses, and counting them would raise the species total to six and
reclassify every group conserved across the four orthopoxviruses as partially
conserved.

The consolidation therefore covers only:

- `cowpox_refseq_Brighton_Red_proteins`
- `mpox_refseq_Rivers_State_proteins`
- `vaccinia_WR_refseq_proteins`
- `varv_refseq_proteins`

## Conservation classes

The classification is derived from the number of species detected, not from a
hard-coded set of four, so adding a species to the input does not require
changing the code:

- `Conserved`: present in all detected species
- `Partially_Conserved`: present in 2 to `n - 1` species
- `Species-Specific`: present in exactly 1 species

## Usage

```bash
python run_analysis.py
```

## Outputs

`--output` names the consolidated CSV; the other three files are derived from
that name and written next to it.

- `outputs/orthogroups_consolidated_analysis.csv`: one row per taxonomic
  combination rather than one per orthogroup. Each row sums, over the
  orthogroups of that combination, how many orthogroups it holds, how many
  proteins in total, how many count as strict orthologs and how many as strict
  paralogs, how many groups are singletons and how many proteins belong to
  them, and the ortholog and paralog counts per species.
- `outputs/orthogroups_consolidated_analysis_classified_protein_ids.txt`: the
  protein IDs grouped by conservation class. This is the file the other steps
  read.
- `outputs/orthogroups_consolidated_analysis_protein_conservation.csv`: one row
  per protein with its class.
- `outputs/orthogroups_consolidated_analysis_expanded_pairwise_conservation.csv`:
  one row per protein and per compared virus.
