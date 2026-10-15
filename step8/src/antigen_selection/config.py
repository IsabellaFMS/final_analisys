from __future__ import annotations

from pathlib import Path

from .models import StepPaths

# BepiPred-3.0 score cutoffs. The high-stringency one is the selection criterion;
# the intermediate one is only used to report the epitopes of an orthologue that
# replaces a truncated member.
HIGH_STRINGENCY = 0.2512
INTERMEDIATE_STRINGENCY = 0.2012

# Step 5 strategic classifications accepted as surface accessible. Group 1 is the
# agreement of SignalP6 and DeepTMHMM; groups 2 and 3 are single-predictor
# positives and are not treated as a consensus.
ACCESSIBLE_GROUPS = ("Group 1: High Confidence Consensus (SP+/THMM+)",)

# An orthologue shorter than this fraction of its orthogroup consensus length is
# treated as truncated. The CrmB orthogroup separates at 0.35 against 0.92 for
# its next shortest member, so the threshold is not sensitive within 0.4-0.5.
TRUNCATION_RATIO = 0.5

# Descriptions that mark an incomplete sequence regardless of its length.
PARTIAL_MARKERS = ("partial", "truncated")

# Columns of step1/sources/orthogroups.tsv that are not part of the analysis.
OUTGROUP_COLUMNS = ("Humanalphaherpesvirus1_proteins", "swinepox_proteins")

SPECIES_BY_COLUMN = {
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "vaccinia_WR_refseq_proteins": "VACV",
    "varv_refseq_proteins": "VARV",
}

OPXV_SPECIES = tuple(SPECIES_BY_COLUMN.values())


def default_paths(root: Path | None = None) -> StepPaths:
    root = root or Path(__file__).resolve().parents[2]
    outputs = root / "outputs"
    return StepPaths(
        root=root,
        outputs=outputs,
        tables=outputs / "tables",
        reports=outputs / "reports",
        bcell=outputs / "bcell",
    )


def sibling_step(paths: StepPaths, step: str, *parts: str) -> Path:
    return paths.root.parent.joinpath(step, *parts)
