# Orthopoxvirus epitope selection pipeline

Eight analysis steps used to select conserved, non-host-homologous B-cell and
T-cell epitope candidates from the proteomes of four orthopoxviruses: _Orthopoxvirus cowpox_
(CPXV), _monkeypox_ (MPXV), _vaccinia_ (VACV) and _variola_ (VARV). The final
products are two candidate tables for MHC-I and MHC-II epitopes and one table
of linear B-cell epitopes, each annotated with the conservation class of its
source protein, and the shortlist of antigens that satisfy every selection
rule.

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
the analysis. The protein classification report produced here is read by steps
2, 5, 6, 7 and 8.

"sources/orthogroups.tsv" is the concatenation of OrthoFinder's assigned
orthogroup table with "orthogroups_unassigned_genes.tsv": 339 rows
("OG0000000" to "OG0000338") covering 1,059 proteins, against 211 rows in
OrthoFinder's file of that name. All 128 identifiers of the unassigned table
also occur in the assigned one; "step1/README.md" explains how the overlap is
resolved.

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

Three properties of this step are worth stating, because each one is a question
a reader is likely to have:

- 36 of the 832 proteins (37 to 86 residues) are absent from
  "step2/outputs/tables/step3_final_consolidated_annotations.tsv": neither
  InterProScan nor eggNOG-mapper returned anything for them. That table
  therefore has 796 rows, while the denominator used throughout is 832, which
  counts those 36 as unannotated.
- Pruning redundant parent terms upstream is statistically neutral, because the
  test propagates counts to ancestor terms. Intersecting the pruned sets and
  intersecting ancestor-closed sets both give 744 terms tested and the same FDR
  for the top term.
- The run log reports 41 obsolete GO identifiers that the deposited
  "go-basic.obo" (release 2024-01-17) no longer contains. GOATOOLS skips them
  without error.

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
filtered FASTA that steps 5 to 8 use. It
also puts back the 45 sequence-duplicated accessions that step 3 removed, each
at its original position, and reports which of them are identical to an
accession the filter excluded.

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
exports by rank, maps each peptide back to its protein, counts how many alleles
bind each epitope, and annotates the conservation class. MHC-I uses a rank
cutoff of 0.75; MHC-II is filtered at 5, 7.5 and 10 and the results are
consolidated.

**Step 8 - final antigen selection.** Applies the selection rules to the step
5, 6 and 7 tables: a B-cell epitope at the high-stringency BepiPred cutoff,
conservation across all four species, surface accessibility by the
SignalP6/DeepTMHMM consensus, one representative per orthogroup chosen on
epitope profile and allelic coverage, and the exclusion of orthologues
truncated relative to their orthogroup consensus length. Two antigens survive.
The step writes a verdict per candidate, so the selection can be audited row by
row. See "step8/README.md".

## Running the pipeline

Each step has its own "run_analysis.py". Steps 1 and 3 have no internal
dependency and can run in any order; step 2 needs the step 1 classification
report and "step1/sources/orthogroups.tsv"; steps 5, 6 and 7 need the step 1
report; step 4 needs the step 3 output; step 7 also needs the step 4 output;
and step 8 needs steps 1, 4, 5, 6 and 7.

```bash
cd step1 && python run_analysis.py && cd ..
cd step2 && python run_analysis.py && cd ..
cd step3 && python run_analysis.py && cd ..
cd step4 && python run_analysis.py && cd ..
cd step5 && python run_analysis.py && cd ..
cd step6 && python run_analysis.py && cd ..
cd step7 && python run_analysis.py && cd ..
cd step8 && python run_analysis.py && cd ..
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
| InterProScan | not recorded in the deposit; the annotation run is dated 21-04-2025 |
| GOATOOLS | pinned in "requirements.txt" |
| SignalP6 | not recorded |
| DeepTMHMM | not recorded |
| BepiPred-3.0 | not recorded |
| IEDB MHC-I and MHC-II | not recorded |

Where a version reads "not recorded", it was not captured at the time of the
run and is not recoverable from the deposited exports, which carry no version
field. The human reference used by the step 4 k-mer comparison is likewise not
recorded beyond the deposited FASTA itself.

## Column names across steps

Step 8 joins the step 5, 6 and 7 tables, which name the same thing differently.
The identifier column is "id" in step 5, "Accession" in step 6 and
"Protein_ID" in steps 7 and 8. Separators also differ: steps 5 and 6 write
semicolon-separated CSV, step 7 comma-separated, and steps 1, 2 and 8
tab-separated. Both are kept as deposited, because the published tables are
cited as they stand.

## Scope

The pipeline selects and annotates candidates, up to and including the final
antigen shortlist. It does not perform docking,
structural modelling, population coverage analysis or any experimental
validation, and it does not run the predictors themselves. Rank and score
cutoffs come from the defaults of each predictor unless stated otherwise in the
step README.
