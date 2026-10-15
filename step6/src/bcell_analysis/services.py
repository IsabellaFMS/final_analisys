from __future__ import annotations

from .config import ExtractionConfig
from .loader import build_conservation_map, load_bepipred_data
from .analyzer import extract_epitopes
from .models import ExecutionSummary, StepPaths
from .writer import export_epitopes_and_fasta


def run_bcell_analysis(paths: StepPaths, config: ExtractionConfig) -> ExecutionSummary:
    proteins = load_bepipred_data(paths.input_file, config.separator)
    epitope_map = extract_epitopes(proteins, config)
    conservation = build_conservation_map(paths.conservation_map)
    artifacts = export_epitopes_and_fasta(epitope_map, paths, conservation)

    return ExecutionSummary(
        input_file=paths.input_file,
        output_root=paths.output_root,
        proteins_loaded=len(proteins),
        accessions_with_epitopes=len(epitope_map),
        epitope_count=sum(len(epitopes) for epitopes in epitope_map.values()),
        epitope_table=artifacts.epitope_table,
        fasta_directory=paths.fasta_dir,
    )
