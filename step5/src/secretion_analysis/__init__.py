from __future__ import annotations

from .workflow import (
    run_comparative_workflow,
    run_deeptmhmm_workflow,
    run_full_workflow,
    run_signalp_workflow,
)

__all__ = [
    "run_comparative_workflow",
    "run_deeptmhmm_workflow",
    "run_full_workflow",
    "run_signalp_workflow",
]

__version__ = "1.1.0"
