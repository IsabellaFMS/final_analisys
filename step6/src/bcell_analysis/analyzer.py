from __future__ import annotations

from typing import Dict, List

from .config import ExtractionConfig
from .models import EpitopeRecord, ProteinRecord


def _mark_as_dummy(
    scores: List[float],
    start: int,
    end: int,
    count: int,
    dummy_score: float,
) -> None:
    
    if count > 3:
        for index in range(start, end):
            scores[index] = dummy_score


def _mask_low_score_runs(scores: List[float], cutoff: float, dummy_score: float) -> List[float]:
    masked = list(scores)
    low_run = 0
    for index, score in enumerate(masked):
        if score < cutoff:
            low_run += 1
            continue
        _mark_as_dummy(masked, index - low_run, index, low_run, dummy_score)
        low_run = 0

    _mark_as_dummy(masked, len(masked) - low_run, len(masked), low_run, dummy_score)
    return masked


def _extend_epitope(
    scores: List[float],
    residues: str,
    peak_index: int,
    max_length: int,
    dummy_score: float,
) -> EpitopeRecord:
    epitope = EpitopeRecord(
        positions=[peak_index],
        residues=[residues[peak_index]],
        scores=[scores[peak_index]],
    )

    while True:
        left = epitope.positions[0] - 1
        right = epitope.positions[-1] + 1
        left_increase = False
        right_increase = False

        if (
            left >= 0
            and len(epitope.residues) < max_length
            and scores[left] != dummy_score
        ):
            left_increase = True
            epitope.positions.insert(0, left)
            epitope.residues.insert(0, residues[left])
            epitope.scores.insert(0, scores[left])

        if (
            right < len(scores)
            and len(epitope.residues) < max_length
            and scores[right] != dummy_score
        ):
            right_increase = True
            epitope.positions.append(right)
            epitope.residues.append(residues[right])
            epitope.scores.append(scores[right])

        if left < 0 and right >= len(scores):
            break
        if not left_increase and not right_increase:
            break
        if (
            left >= 0
            and right < len(scores)
            and scores[left] == dummy_score == scores[right]
        ):
            break

    return epitope


def extract_epitopes(
    proteins: List[ProteinRecord], config: ExtractionConfig
) -> Dict[str, List[EpitopeRecord]]:
    
    results: Dict[str, List[EpitopeRecord]] = {}
    dummy_score = config.cutoff - 1e19

    for protein in proteins:
        scores = _mask_low_score_runs(protein.scores, config.cutoff, dummy_score)
        epitopes: List[EpitopeRecord] = []

        while scores and max(scores) >= config.cutoff:
            peak_index = scores.index(max(scores))
            epitope = _extend_epitope(
                scores=scores,
                residues=protein.residues,
                peak_index=peak_index,
                max_length=config.max_length,
                dummy_score=dummy_score,
            )

            if len(epitope.residues) >= config.min_length:
                epitopes.append(epitope)

            for position in epitope.positions:
                scores[position] = dummy_score

        if epitopes:
            results[protein.accession] = epitopes

    return results
