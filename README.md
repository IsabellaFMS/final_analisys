# Orthopoxvirus epitope selection pipeline

Seven analysis steps used to select conserved, non-host-homologous B-cell and
T-cell epitope candidates from the proteomes of four orthopoxviruses: cowpox
(CPXV), monkeypox (MPXV), vaccinia (VACV) and variola (VARV). The final
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
Swiss-Prot entries was run outside this repository; the step consumes the
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

## Data flow

| Step | Main input | Comes from | Main output |
|------|-----------|------------|-------------|
| 1 | `orthogroups.tsv`, `orthogroups_unassigned_genes.tsv` | OrthoFinder | `orthogroups_consolidated_analysis_classified_protein_ids.txt` |
| 2 | eggNOG and InterPro exports, step 1 classification | external tools, step 1 | enrichment tables and GO matrices |
| 3 | `human.fasta`, `viral_proteins.fasta` | Swiss-Prot, RefSeq | `clean_human_proteome.fasta`, `clean_viral_proteins.fasta`, `removed_viral_proteins.fasta` |
| 4 | `viral_hits.tsv`, `clean_viral_proteins.fasta`, the step 3 removal report | external k-mer run, step 3 | `safe_sequences.fasta`, `exclusion_ids.txt` |
| 5 | SignalP6 and DeepTMHMM exports, step 1 classification | external tools, step 1 | secretion and topology tables |
| 6 | BepiPred-3.0 export, step 1 classification | external tool, step 1 | `epitopes.csv` |
| 7 | IEDB exports, `safe_sequences.fasta`, step 1 classification | external tool, steps 1 and 4 | `mhci_final.csv`, `mhcii_final.csv` |

Step 1 is the only producer of the conservation report. Steps 5, 6 and 7 read
it directly from "step1/outputs/" and fall back to a local copy under
"sources/" when step 1 has not been run. Step 2 reads its own copy under the
shorter name "classified_protein_ids.txt", kept from an earlier version of the
analysis.

Files under any "sources/" directory keep the names they were given by the run
that produced them, which is why some of them are in Portuguese and why the
same classification report appears under two different names.

Step 3 was introduced after the pipeline had already been run once. Removing the
45 duplicated viral sequences upstream would take "safe_sequences.fasta" from
191 to 146 records and invalidate the SignalP6, DeepTMHMM, BepiPred-3.0 and
IEDB exports stored under "sources/" in steps 5, 6 and 7, which were produced
on the 191-record proteome. Step 4 therefore puts the 45 back, each at the
position it held in the raw file, and writes a FASTA identical in content and
in order to the one those predictors were run on. Deduplication costs no
epitope content, the sequences being identical, and it changes no conservation
annotation: all 45 removed accessions carry the same class as the accession
kept in their place. The four predictors do not need to be run again.

One caveat travels with that decision, and the step 4 execution summary counts
it: **33 of the 45 are identical to an accession that step 4 excluded for
human homology**, and are back in the filtered proteome despite sharing the
human 12-mers that got their twin excluded. They were absent from
"viral_hits.tsv" only because the external k-mer comparison collapsed repeated
sequences and reported each hit through a single representative. This is a
property the 191-record proteome always had — deduplication had hidden it, not
created it — and the 33 are flagged rather than silently included, so that any
candidate derived from them can be traced. The other 12 pass the filter on
their own.

## External tools

| Tool | Version | Used by |
|------|---------|---------|
| OrthoFinder | not recorded | step 1 |
| eggNOG-mapper | v2.1.12, DIAMOND mode, restricted to Viruses, E-value <= 1e-5, identity >= 40%, coverage >= 20% | step 2 |
| InterProScan | v5.75-106.019 | step 2 |
| GOATOOLS | not recorded | step 2 |
| R | v4.3.3 | step 2 |
| SignalP | 6 | step 5 |
| DeepTMHMM | not recorded | step 5 |
| BepiPred | 3.0 | step 6 |
| NetMHCpan / NetMHCIIpan | via IEDB, "netmhciipan_el" method for MHC-II | step 7 |

Versions marked as not recorded were not written down at the time of the run.

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
