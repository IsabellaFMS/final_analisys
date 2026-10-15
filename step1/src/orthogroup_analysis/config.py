from __future__ import annotations

import re

from orthogroup_analysis.models import SpeciesConfig


ORTHOGROUP_COLUMN = "Orthogroup"

# Human alphaherpesvirus 1 and swinepox are present in the OrthoFinder output
# but are not orthopoxviruses, so they take no part in the conservation classes.
EXCLUDED_HEADERS = {
    "Humanalphaherpesvirus1_proteins",
    "swinepox_proteins",
}

ALIAS_OVERRIDES = {
    "Humanalphaherpesvirus1_proteins": "HSV1",
    "cowpox_refseq_Brighton_Red_proteins": "CPXV",
    "mpox_refseq_Rivers_State_proteins": "MPXV",
    "swinepox_proteins": "SWPV",
    "vaccinia_WR_refseq_proteins": "VACV",
    "varv_refseq_proteins": "VARV",
}

DISPLAY_NAME_OVERRIDES = {
    "Humanalphaherpesvirus1_proteins": "Human alphaherpesvirus 1",
    "cowpox_refseq_Brighton_Red_proteins": "Cowpox virus Brighton Red",
    "mpox_refseq_Rivers_State_proteins": "Mpox virus Rivers State",
    "swinepox_proteins": "Swinepox virus",
    "vaccinia_WR_refseq_proteins": "Vaccinia virus WR",
    "varv_refseq_proteins": "Variola virus",
}


def build_species_config(headers: list[str]) -> list[SpeciesConfig]:
    species_headers = [
        header
        for header in headers
        if header != ORTHOGROUP_COLUMN and header not in EXCLUDED_HEADERS
    ]

    return [
        SpeciesConfig(
            header=header,
            alias=ALIAS_OVERRIDES.get(header, build_alias(header)),
            display_name=DISPLAY_NAME_OVERRIDES.get(header, build_display_name(header)),
        )
        for header in species_headers
    ]


def build_alias(header: str) -> str:
    cleaned = re.sub(r"_proteins$", "", header)
    parts = [part for part in re.split(r"[_\W]+", cleaned) if part]
    if not parts:
        return header.upper()

    # Initials of the first six parts
    capitals = "".join(part[0].upper() for part in parts[:6])
    return capitals or header.upper()


def build_display_name(header: str) -> str:
    cleaned = re.sub(r"_proteins$", "", header)
    cleaned = cleaned.replace("_", " ").strip()
    return cleaned or header
