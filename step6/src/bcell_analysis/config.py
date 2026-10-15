from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ExtractionConfig:
    cutoff: float
    min_length: int = 20
    max_length: int = 30
    separator: str = ","

    def __post_init__(self) -> None:
        if self.min_length <= 0:
            raise ValueError("min_length must be positive.")
        if self.max_length < self.min_length:
            raise ValueError("max_length must be greater than or equal to min_length.")
        if not self.separator:
            raise ValueError("separator must not be empty.")
