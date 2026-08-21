# Orthopoxvirus epitope selection pipeline

Seven analysis steps used to select conserved, non-host-homologous B-cell and
T-cell epitope candidates from the proteomes of four orthopoxviruses: _Orthopoxvirus cowpox_
(CPXV), _monkeypox_ (MPXV), _vaccinia_ (VACV) and _variola_ (VARV). The final
products are two candidate tables for MHC-I and MHC-II epitopes and one table
of linear B-cell epitopes, each annotated with the conservation class of its
source protein.

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
2, 5, 6 and 7.

**Step 2 - functional annotation and GO enrichment.** Merges the eggNOG-mapper
and InterProScan annotations at protein level, removes redundant parent GO
terms when a more specific descendant is present, and runs GO enrichment per
conservation group with GOATOOLS, using Fisher's exact test with
Benjamini-Hochberg correction at FDR 0.05. The proteome of each virus is the
background population. An R script builds the interactive heatmap of the top 10
GO terms per category, with the four proteomes side by side.

**Step 3 - proteome cleaning.** Drops duplicated sequences and sequences
shorter than the k-mer size from the raw human and viral proteomes, since
neither can produce a usable comparison. Duplicates are resolved by sequence and
the first occurrence in file order is kept. Everything removed is logged with
its reason. The two cleaned FASTAs are what the k-mer comparison consumes, and
the viral one is the input of step 4.

**Step 4 - host homology filtering.** Removes viral proteins that share
identical 12-residue stretches with the human proteome, to limit
cross-reactivity in downstream candidates. The k-mer comparison against human
UniProt entries was run outside this repository; the step consumes the
resulting hit table and writes the filtered FASTA that steps 5 to 7 use. It
also puts back the 45 sequence-duplicated accessions that step 3 removed, each
at its original position, and reports which of them are identical to an
accession the filter excluded.

**Step 5 - secretion and topology.** Combines the SignalP6 signal-peptide
predictions with the DeepTMHMM transmembrane topology and annotates each
protein with its conservation class, separating extracellularly accessible
proteins from intracellular ones.

**Step 6 - linear B-cell epitopes.** Reads the residue-level BepiPred-3.0
export and extracts contiguous regions above the score cutoff, keeping those
between 20 and 30 residues.

**Step 7 - T-cell epitopes.** Filters the IEDB MHC-I and MHC-II prediction
exports by rank, maps each peptide back to its protein, counts how many alleles
bind each epitope, and annotates the conservation class. MHC-I uses a rank
cutoff of 0.75; MHC-II is filtered at 5, 7.5 and 10 and the results are
consolidated.

## Running the pipeline

Each step has its own "run_analysis.py". Steps 1 and 4 have no internal
dependency and can run in any order; steps 2, 5, 6 and 7 need the step 1
classification report, and step 7 also needs the step 4 output.

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
"dplyr" and "tidyr".

See the "README.md" inside each step for its parameters and outputs.

## Scope

The pipeline selects and annotates candidates. It does not perform docking,
structural modelling, population coverage analysis or any experimental
validation, and it does not run the predictors themselves. Rank and score
cutoffs come from the defaults of each predictor unless stated otherwise in the
step README.
