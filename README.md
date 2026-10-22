# Orthopoxvirus epitope selection pipeline

Seven analysis steps used to select conserved, non-host-homologous B-cell and
T-cell epitope candidates from the proteomes of four orthopoxviruses: _Orthopoxvirus cowpox_
(CPXV), _monkeypox_ (MPXV), _vaccinia_ (VACV) and _variola_ (VARV).

Each step is a self-contained project with its own "sources/", "outputs/" and
"run_analysis.py". The predictions themselves are not produced here: the
external tools listed below were run separately, and their exports are the
inputs of the corresponding step.

## Steps

**Step 1 - conservation map.** Reads the OrthoFinder orthogroup tables and
classifies each orthogroup by how many of the four species it contains:
conserved (all four), partially conserved (two or three), species-specific
(one). Two species columns present in the OrthoFinder output,
"Humanalphaherpesvirus1_proteins" and "swinepox_proteins", are excluded from
the analysis.

"sources/orthogroups.tsv" is the concatenation of OrthoFinder's assigned
orthogroup table with "orthogroups_unassigned_genes.tsv": 339 rows
("OG0000000" to "OG0000338") covering 1,059 proteins, against 211 rows in
OrthoFinder's file of that name.

**Step 2 - functional annotation and GO term distribution.** Merges the
eggNOG-mapper and InterProScan annotations at protein level, removes redundant
parent GO terms when a more specific descendant is present, and tests each
conservation group with GOATOOLS, using Fisher's exact test with
Benjamini-Hochberg correction at FDR 0.05. The background population is the set
of 483 proteins carrying at least one GO term, pooled across the four
proteomes. The test runs at two levels of analysis: the protein, which is
exploratory because orthologues of one family are not independent observations,
and the orthogroup, which is the level reported in the manuscript. No term is
significant at FDR 0.05 at the orthogroup level. An R script builds the
interactive heatmap of the top 10 GO terms per category, with the four
proteomes side by side.

- 36 of the 832 proteins (37 to 86 residues) are absent from
  "step2/outputs/tables/final_consolidated_annotations.tsv": neither
  InterProScan nor eggNOG-mapper returned anything for them. That table
  therefore has 796 rows, while the denominator used throughout is 832, which
  counts those 36 as unannotated.
- Pruning redundant parent terms upstream is statistically neutral, because the
  test propagates counts to ancestor terms. Intersecting the pruned sets and
  intersecting ancestor-closed sets both give 744 terms tested and the same FDR
  for the top term.
- 45 of the GO identifiers in the annotations are flagged "is_obsolete" in the
  deposited "go-basic.obo" (release 2024-01-17). The file still defines them,
  so they are tested as deposited and GOATOOLS raises no error.

**Step 3 - proteome cleaning.** Drops duplicated sequences and sequences
shorter than the k-mer size from the raw human and viral proteomes, since
neither can produce a usable comparison. Duplicates are resolved by sequence and
the first occurrence in file order is kept. Everything removed is logged with
its reason. The two cleaned FASTAs are what the k-mer comparison consumes, and
the viral one is the input of step 4.

**Step 4 - host homology filtering.** Removes viral proteins that contain a
12-residue window with 9 or more residues identical to a human 12-mer, to limit
cross-reactivity in downstream candidates. Of the 19,709 hit rows, 1,078 are
exact 12/12 matches and 13,995 match at 9/12, over 641 distinct viral
proteins. The k-mer comparison against human UniProt entries was run outside
this repository, with a window of 12 residues and an acceptance threshold of 9
identical residues; the step consumes the resulting hit table and writes the
filtered FASTA that steps 5 to 7 use. It
also puts back the 45 sequence-duplicated accessions that step 3 removed, each
at its original position.

**Step 5 - secretion and topology.** Combines the SignalP6 signal-peptide
predictions with the DeepTMHMM transmembrane topology and annotates each
protein with its conservation class, separating extracellularly accessible
proteins from intracellular ones.

**Step 6 - linear B-cell epitopes.** Reads the residue-level BepiPred-3.0
export, masks runs of more than three consecutive sub-threshold residues,
anchors each candidate at the highest remaining score and grows it
symmetrically, keeping regions of 20 to 30 residues. "step6/README.md"
describes the procedure in full.

**Step 7 - T-cell epitopes.** Filters the IEDB MHC-I and MHC-II prediction
exports by rank, maps each epitope back to its protein, counts how many alleles
bind each epitope, and annotates the conservation class. MHC-I uses a rank
cutoff of 0.75. MHC-II is filtered at rank 10: the code applies 5, 7.5 and 10
and unions the three, but the first two sets are contained in the third, so the
consolidated set is the rank 10 set and the three cutoffs have no effect on it.

## Running the pipeline

Each step has its own "run_analysis.py". Steps 1 and 3 have no internal
dependency and can run in any order; step 2 needs
"step1/sources/orthogroups.tsv", which it reads to map each protein to its
orthogroup; steps 5, 6 and 7 need the step 1 classification report.

```bash
cd step1 && python run_analysis.py && cd ..
cd step2 && python run_analysis.py && cd ..
cd step3 && python run_analysis.py && cd ..
cd step4 && python run_analysis.py && cd ..
cd step5 && python run_analysis.py && cd ..
cd step6 && python run_analysis.py && cd ..
cd step7 && python run_analysis.py && cd ..
```

Python 3.11 or newer, with the packages in "requirements.txt". The R step needs
"dplyr", "tidyr", "plotly" and "htmlwidgets", and "Rscript" on PATH or
"OPXV_RSCRIPT" pointing at it.

See the "README.md" inside each step for its parameters and outputs.

## Tool versions

The predictions are not produced here. The external tools were run separately,
in these versions:

| Tool | Version |
|---|---|
| OrthoFinder | v2.5.5, MCL inflation 1.7 |
| eggNOG-mapper | v2.1.12 |
| InterProScan | v5.73-104.0 |
| GOATOOLS | v1.6.4 |
| SignalP6 | v6.0 |
| DeepTMHMM | v1.0 |
| BepiPred-3.0 | v3.0 |
| IEDB MHC-I and MHC-II | v4.1 and v4.0 |
