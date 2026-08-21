from __future__ import annotations

from .models import ExecutionSummary


def format_execution_summary(summary: ExecutionSummary) -> str:
    return "\n".join(
        [
            "Execution summary",
            f"Input file: {summary.input_file}",
            f"Output root: {summary.output_root}",
            f"Proteins loaded: {summary.proteins_loaded}",
            f"Accessions with epitopes: {summary.accessions_with_epitopes}",
            f"Total epitopes: {summary.epitope_count}",
            f"Epitope table: {summary.epitope_table}",
            f"FASTA directory: {summary.fasta_directory}",
        ]
    )
