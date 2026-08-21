from __future__ import annotations

import argparse
from pathlib import Path

import sys

SRC_DIR = Path(__file__).resolve().parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from orthogroup_analysis.analyzer import analyze_orthogroups, build_classified_protein_groups, build_expanded_conservation_rows, build_protein_conservation_rows
from orthogroup_analysis.loader import load_records
from orthogroup_analysis.writer import write_classified_protein_ids_txt, write_consolidated_csv, write_expanded_conservation_csv, write_protein_conservation_csv


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUTS = [
    PROJECT_ROOT / "sources" / "orthogroups.tsv",
    PROJECT_ROOT / "sources" / "orthogroups_unassigned_genes.tsv",
]
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "orthogroups_consolidated_analysis.csv"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=("Consolidates orthogroup analysis from two or more TSV files.")
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        help=(
            "Input TSV files. Defaults to sources/orthogroups.tsv and "
            "sources/orthogroups_unassigned_genes.tsv."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output CSV file. Defaults to outputs/orthogroups_consolidated_analysis.csv.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_paths = args.inputs or DEFAULT_INPUTS
    output_path = args.output or DEFAULT_OUTPUT

    species_config, records = load_records(input_paths)
    consolidated_rows = analyze_orthogroups(records, species_config)
    expanded_rows = build_expanded_conservation_rows(records, species_config)
    protein_conservation_rows = build_protein_conservation_rows(records, species_config)
    classified_protein_groups = build_classified_protein_groups(records, species_config)

    expanded_output_path = build_expanded_output_path(output_path)
    protein_conservation_output_path = build_protein_conservation_output_path(output_path)
    classified_protein_ids_output_path = build_classified_protein_ids_output_path(output_path)
    write_consolidated_csv(output_path, consolidated_rows, species_config)
    write_expanded_conservation_csv(expanded_output_path, expanded_rows)
    write_protein_conservation_csv(
        protein_conservation_output_path,
        protein_conservation_rows,
    )
    write_classified_protein_ids_txt(
        classified_protein_ids_output_path,
        classified_protein_groups,
    )

    print(f"Files read: {', '.join(path.name for path in input_paths)}")
    print(f"Orthogroups processed: {len(records)}")
    print(f"Consolidated taxonomic combinations: {len(consolidated_rows)}")
    print(f"Consolidated output written to: {output_path}")
    print(f"Expanded pairwise output written to: {expanded_output_path}")
    print(f"Protein conservation output written to: {protein_conservation_output_path}")
    print(f"Classified protein IDs written to: {classified_protein_ids_output_path}")


def build_expanded_output_path(output_path: Path) -> Path:
    suffix = output_path.suffix or ".csv"
    return output_path.with_name(
        f"{output_path.stem}_expanded_pairwise_conservation{suffix}"
    )


def build_protein_conservation_output_path(output_path: Path) -> Path:
    suffix = output_path.suffix or ".csv"
    return output_path.with_name(
        f"{output_path.stem}_protein_conservation{suffix}"
    )


def build_classified_protein_ids_output_path(output_path: Path) -> Path:
    
    return output_path.with_name(
        f"{output_path.stem}_classified_protein_ids.txt"
    )


if __name__ == "__main__":
    main()
