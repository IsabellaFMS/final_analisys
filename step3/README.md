# Step 3: proteome cleaning

Prepares the raw human and viral proteomes for the k-mer comparison and for the
filtering of step 4, by dropping the records that cannot take part in that
comparison: sequences that are duplicated, which would only multiply the same
k-mers and the same result rows, and sequences shorter than the k-mer size,
which produce no k-mer at all.

Everything that is removed is logged with its reason, so the cleaning can be
audited.

## Inputs

In `sources/`:

- `human.fasta`, the human proteome as downloaded from UniProt.
- `viral_proteins.fasta`, the RefSeq proteomes of the four orthopoxviruses.

Both are the raw downloads. They keep the names they were given by the run that
produced them, which is why one of them is in Portuguese.

## How cleaning works

Each proteome is treated independently, in this order:

1. **Too short.** A sequence shorter than `--kmer-size` is dropped, with reason
   `Short`. The default is 12, the k used in the comparison downstream. A
   sequence of exactly 12 residues yields one 12-mer and is kept.
2. **Duplicated.** Among the remaining records, sequences are grouped by exact
   match and only the **first occurrence in file order** is kept. The others are
   dropped with reason `Duplicate` and a `kept_instead` pointing at the survivor.

Short records are screened out first, so one can never be picked as the survivor
of a group of identical sequences.

The order of the surviving records is the order of the input file, which matters
for any downstream tool that indexes by position.

Deduplication is by **sequence**, not by header. In the current viral proteome
that collapses 42 groups, and **19 of them span species** the same protein in,
say, vaccinia and cowpox with an identical sequence. The accession that survives
is the one from whichever species comes first in the file, and the ortholog of
the other species leaves the analysis. No epitope content is lost, since the
sequences are identical, but the conservation annotation of the dropped
accession is.

Records that share an accession but not a sequence are reported as a warning and
nothing is removed: step 4 matches by accession, and a repeated one would make
that match ambiguous.

## Usage

```bash
python run_analysis.py
```

## Outputs

In `outputs/`:

- `fasta/clean_human_proteome.fasta` and `fasta/clean_viral_proteins.fasta`,
  the two cleaned proteomes.
- `reports/removed_human_proteins.csv` and `reports/removed_viral_proteins.csv`,
  one row per removed record with its reason, the accession kept in its place
  when the reason is `Duplicate`, and `source_index`, the position the record
  held in the raw file.
- `fasta/removed_viral_proteins.fasta`, the viral records that were removed.
- `reports/cleaning_summary.txt`, the run counts.

The removed-viral FASTA and the `source_index` column exist for step 4, which
reads both to put the 45 duplicated viral accessions back into the filtered
proteome, each at its original position. See the step 4 README.

## Notes

- The k-mer comparison itself is not run here. This step produces the two FASTAs
  that it consumes; its output table is the input of step 4.

