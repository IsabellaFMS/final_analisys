# Step 7: T-cell epitope screening

Filters the IEDB MHC-I and MHC-II prediction exports by rank, maps each epitope
back to its source protein, counts how many alleles bind each epitope and
annotates the conservation class from step 1. The proteome screened is the one
step 4 produced.

## Inputs

Read from `sources/`:

- `mhci/full_predictions/prediction_*.csv`
- `mhcii/DRB1/DRB1_predict_result.csv` and
  `mhcii/DQA1/DQA1_predict_result.csv`
- `safe_sequences.fasta`, the filtered proteome from step 4. Epitopes are tied
  back to proteins by the position of the sequence in this file, through the
  `seq_num` column, so the FASTA has to be the same one submitted to IEDB.
- `orthogroups_consolidated_analysis_classified_protein_ids.txt`, the step 1
  report


## Usage

```bash
python run_analysis.py
```

## Outputs

`outputs/mhci_final.csv` and `outputs/mhcii_final.csv`, one row per protein,
allele and peptide, with `Allele Count` giving the number of distinct alleles
that bind the peptide and `Alleles` listing them.

MHC-I is filtered at percentile rank 0.75. MHC-II is filtered at rank 10: the
code applies 5, 7.5 and 10 and unions the three, but the first two sets are
contained in the third, so the consolidated set is the rank 10 set.
