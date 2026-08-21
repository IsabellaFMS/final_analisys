# Step 7: T-cell epitope screening

Filters the IEDB MHC-I and MHC-II prediction exports by rank, maps each epitope
back to its source protein, counts how many alleles bind each epitope and
annotates the conservation class from step 1. The proteome screened is the one
step 4 produced.

## Inputs

Read from "sources/":

- "mhci/full_predictions/prediction_*.csv", the MHC-I predictions.
- "mhcii/DRB1/DRB1_predict_result.csv" and
  "mhcii/DQA1/DQA1_predict_result.csv", the MHC-II predictions. The
  parameter files next to them record the alleles and the IEDB method used
  ("netmhciipan_el").
- "safe_sequences.fasta", the filtered proteome from step 4. Epitopes are tied
  back to proteins by the position of the sequence in this file, through the
  "seq_num" column, so the FASTA has to be the same one submitted to IEDB.
- "orthogroups_consolidated_analysis_classified_protein_ids.txt", the step 1
  report. It is read from "step1/outputs/" when available and from this local
  copy otherwise.


## Usage

```bash
python run_analysis.py
```

## Outputs

"outputs/mhci_final.csv" and "outputs/mhcii_final.csv".