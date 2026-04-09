# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Public exports for the FCC protein-folding package.

The package exposes qiskit-dependent classes lazily so that lightweight
submodules such as :mod:`fcc.classical_utils` can still be imported in
environments where qiskit is not installed.
"""

from importlib import import_module


_EXPORT_MAP = {
    "ProteinFoldingProblem": ("fcc_protein_folding_problem", "ProteinFoldingProblem"),
    "ProteinFoldingResult": ("fcc_protein_folding_result", "ProteinFoldingResult"),
    "Peptide": ("fcc_peptide", "Peptide"),
    "PenaltyParameters": ("fcc_penalty_parameters", "PenaltyParameters"),
    "HPInteraction": ("fcc_hp_interaction", "HPInteraction"),
    "MiyazawaJerniganInteraction": ("fcc_mj_interaction", "MiyazawaJerniganInteraction"),
    "ProteinSolver": ("fcc_solver", "ProteinSolver"),
    "ProteinShapeDecoder": ("fcc_protein_shape", "ProteinShapeDecoder"),
    "ProteinShapeFileGen": ("fcc_protein_shape", "ProteinShapeFileGen"),
}

__all__ = list(_EXPORT_MAP)


def __getattr__(name: str):
    if name not in _EXPORT_MAP:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORT_MAP[name]
    module = import_module(f".{module_name}", __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(list(globals().keys()) + __all__)
