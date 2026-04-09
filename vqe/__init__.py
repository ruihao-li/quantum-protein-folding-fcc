# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Public exports for the VQEC optimization helpers.

The VQE helpers depend on Qiskit. Import them lazily so importing :mod:`vqe`
does not immediately fail in environments that only need package metadata or
discovery.
"""

from importlib import import_module


_EXPORT_MAP = {
    "PerturbedPrimalDualOpt": ("vqec_optimization", "PerturbedPrimalDualOpt"),
    "OptimisticGDAOpt": ("vqec_optimization", "OptimisticGDAOpt"),
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
