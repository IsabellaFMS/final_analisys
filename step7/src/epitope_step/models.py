from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class WorkflowResult:
    workflow: str
    artifacts: list[Path] = field(default_factory=list)
    metrics: dict[str, int] = field(default_factory=dict)
    final_frame: pd.DataFrame | None = None
