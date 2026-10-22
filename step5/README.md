# Step 5: secretion and topology

Crosses the SignalP6 signal peptide predictions with the DeepTMHMM
topology and annotates the result with the conservation label
from step 1. The point is to separate the proteins that are accessible outside
the cell from the ones that are not.

## Inputs

- `sources/signalp6/prediction_results.txt` and
  `sources/signalp6/region_output.gff3`
- `sources/deeptmhmm/TMRs.gff3`
- `sources/orthogroups_consolidated_analysis_classified_protein_ids.txt`

Both predictors were run outside this repository, on the proteome produced by
step 4.

## Usage

```bash
python run_analysis.py 
```

## Outputs

- `outputs/tables/signalp6_results_final.csv`
- `outputs/tables/deeptmhmm_summary_final.csv`
- `outputs/tables/comparative_results_signalp6_deep_final.csv`


