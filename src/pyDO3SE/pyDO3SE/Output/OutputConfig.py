from dataclasses import dataclass, field
from typing import List

from .Output_Shape import default_output_fields

@dataclass
class OutputConfig:

    fields: List[str] = field(default_factory=lambda: default_output_fields)
    log_multilayer: bool = False
