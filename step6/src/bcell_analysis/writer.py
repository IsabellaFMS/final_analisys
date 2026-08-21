from __future__ import annotations

import csv
import statistics
from pathlib import Path
from typing import Dict, List

from .loader import UNCLASSIFIED_LABEL, protein_accession
from .models import EpitopeRecord, OutputArtifacts, StepPaths


def _sanitize_filename(name: str) -> str:
    sanitized = name.replace(" ", "_")
    for symbol in '<>:"/\\|?*[]{}()&%$#@!;':
        sanitized = sanitized.replace(symbol, "")
    while "__" in sanitized:
        sanitized = sanitized.replace("__", "_")
    return sanitized.strip(" ._")


def export_epitopes_and_fasta(
    epitope_map: Dict[str, List[EpitopeRecord]],
    paths: StepPaths,
    conservation: Dict[str, str],
) -> OutputArtifacts:
    
    paths.tables_dir.mkdir(parents=True, exist_ok=True)
    paths.fasta_dir.mkdir(parents=True, exist_ok=True)

    epitope_table = paths.tables_dir / "epitopes.csv"
    fasta_files: List[Path] = []

    with epitope_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(
            [
                "Accession",
                "ID",
                "Epitope",
                "Average Score",
                "Max Score",
                "Median Score",
                "Start position",
                "End position",
                "Length",
                "Classification (conservation)",
            ]
        )

        for accession, epitopes in epitope_map.items():
            classification = conservation.get(protein_accession(accession), UNCLASSIFIED_LABEL)
            fasta_path = paths.fasta_dir / f"{_sanitize_filename(accession)}.fasta"
            fasta_files.append(fasta_path)

            with fasta_path.open("w", encoding="utf-8", newline="") as fasta_handle:
                for index, epitope in enumerate(epitopes, start=1):
                    writer.writerow(
                        [
                            accession,
                            index,
                            epitope.residue_sequence,
                            statistics.mean(epitope.scores),
                            max(epitope.scores),
                            statistics.median(epitope.scores),
                            epitope.start_position,
                            epitope.end_position,
                            epitope.length,
                            classification,
                        ]
                    )
                    fasta_handle.write(
                        f">{accession} | epitope {index} | start {epitope.start_position} | "
                        f"end {epitope.end_position}\n"
                    )
                    fasta_handle.write(f"{epitope.residue_sequence}\n")

    return OutputArtifacts(epitope_table=epitope_table, fasta_files=fasta_files)
