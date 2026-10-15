# Step 6: linear B-cell epitopes

Extracts linear B-cell epitope candidates from the residue-level BepiPred-3.0
export, and labels each one with the conservation class of its protein.

## Inputs

- "sources/0.1512_raw_output.csv", the BepiPred-3.0 export. One row per
  residue, with the accession, the residue, the BepiPred-3.0 score and the
  linear epitope score. The "0.1512" in the name is the threshold that was set
  when the export was generated.
- "sources/orthogroups_consolidated_analysis_classified_protein_ids.txt", the
  step 1 report. This is the fallback copy: the resolver prefers
  "../step1/outputs/", and "--classification-file" overrides both. The report
  is required, and the step stops with an error if it is missing rather than
  writing an unlabelled table.


## How extraction works

1. Runs of more than three consecutive residues below the cutoff are blanked
   out. Shorter dips are left alone, so a single low-scoring residue does not
   split a reactive region in two.
2. The highest remaining score is taken as a seed and the region grows on both
   sides until it hits a blanked position on both sides or reaches
   "--length-max".
3. The region is kept if it is at least "--length-min" long, then its positions
   are blanked and the search repeats on what is left.

Because growth is not conditioned on the score, a candidate can include
residues below the cutoff, provided the dip was short enough to survive step 1.

## Usage

```bash
python run_analysis.py
```

## Outputs

- "outputs/tables/epitopes.csv", one row per epitope: "Accession", "ID",
  "Epitope", the average, maximum and median score, "Start position",
  "End position", "Length" and "Classification (conservation)", the class the
  step 1 report gives the source protein. The column carries the same name it
  has in the step 7 tables, and a protein absent from the report would read
  "Unclassified" — none does in the current run.
- "outputs/fasta/*.fasta", one file per protein with candidates.
