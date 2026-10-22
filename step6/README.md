# Step 6: linear B-cell epitopes

Extracts linear B-cell epitope candidates from the residue-level BepiPred-3.0
export, and labels each one with the conservation class of its protein.

## Inputs

- `sources/0.1512_raw_output.csv`, the BepiPred-3.0 export. One row per
  residue, with the accession, the residue, the BepiPred-3.0 score and the
  linear epitope score.
- `sources/orthogroups_consolidated_analysis_classified_protein_ids.txt`


## How extraction works

1. Runs of more than three consecutive residues below the cutoff are blanked
   out. Shorter dips are left alone, so a single low-scoring residue does not
   split a reactive region in two.
2. The highest remaining score is taken as a seed and the region grows on both
   sides until it hits a blanked position on both sides or reaches
   `--length-max`.
3. The region is kept if it is at least `--length-min` long, then its positions
   are blanked and the search repeats on what is left.

Because growth is not conditioned on the score, a candidate can include
residues below the cutoff.

## Usage

```bash
python run_analysis.py
```

## Outputs

- `outputs/tables/epitopes.csv`
- `outputs/fasta/*.fasta`
- `outputs/tables/epitopes_0.2012.csv` and `outputs/tables/epitopes_0.2512.csv`,
  the same table at the two higher cutoffs, written by rerunning with
  `--cutoff`. `epitopes.csv` is the default 0.1512 run and holds 473 epitopes;
  the other two hold 168 and 40. The 0.2512 table is the one the manuscript
  reports.

```bash
python run_analysis.py --cutoff 0.2512 --output-dir <dir>
```

The output directory is not cleared between runs, so a rerun at a different
cutoff into the same directory leaves the previous run's FASTA files behind.
Pass a fresh `--output-dir` for each cutoff.
