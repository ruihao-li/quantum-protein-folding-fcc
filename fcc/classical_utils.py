# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Utility functions relevant for classical search methods."""

import os
import sys
from pathlib import Path


def _relabel_turns(turns: str) -> str:
    """Relabels the turns from the classical search by Frank to the ones used in
    the quantum encoding by Ruihao. Note that the turns from the classical
    search are read from right to left. After relabeling, the turns will be read
    from left to right."""
    mapping = {
        "0": "0",
        "1": "2",
        "2": "4",
        "3": "6",
        "4": "8",
        "5": "a",
        "6": "1",
        "7": "3",
        "8": "5",
        "9": "7",
        "a": "9",
        "b": "b",
    }
    try:
        new_turns = "".join(mapping[turn.lower()] for turn in turns)
    except KeyError as exc:
        raise ValueError(
            f"Unsupported classical-search turn symbol: {exc.args[0]!r}."
        ) from exc
    return new_turns[::-1]  # Reverse the string to read from left to right


def _construct_results_file_path(file_name: str, results_dir: str | None = None) -> str:
    """Construct an absolute path to a classical-search results file.

    The repository layout stores the result files under
    ``<repo_root>/classical_search/results``, while package installs place them
    under ``<sys.prefix>/share/quantum-protein-folding-fcc/classical_search/results``.
    """
    if results_dir is not None:
        return os.path.normpath(os.path.join(results_dir, file_name))

    repo_root = Path(__file__).resolve().parent.parent
    candidate_dirs = [
        repo_root / "classical_search" / "results",
        Path(sys.prefix)
        / "share"
        / "quantum-protein-folding-fcc"
        / "classical_search"
        / "results",
    ]

    for candidate_dir in candidate_dirs:
        candidate_file = candidate_dir / file_name
        if candidate_file.exists():
            return os.path.normpath(str(candidate_file))

    searched_dirs = ", ".join(str(path) for path in candidate_dirs)
    raise FileNotFoundError(
        f"Could not find classical-search results file '{file_name}'. Searched: {searched_dirs}"
    )


def load_top_cls_solns(
    file_name: str, results_dir: str | None = None
) -> list[tuple[str, float]]:
    """
    Loads the top solutions from the classical search results file
    (`topobj_<sequence>.txt`).

    Args:
        file_name (str): The name of the file containing the classical search results.
        results_dir (str | None): Optional directory containing classical search
            result files. If omitted, defaults to
            `<repo_root>/classical_search/results`.

    Returns:
        list[tuple[str, float]]: A list of tuples where each tuple contains a
        string representing the turns and a float representing the energy of
        that configuration.
    """
    top_cls_turns: list[tuple[str, float]] = []
    file_path = _construct_results_file_path(file_name, results_dir)
    with open(file_path, "r") as file:
        data = file.read()
        # Split by line
        data_by_line = data.split("\n")
        for line in data_by_line:
            if len(line.split("\t")) == 2:
                [energy, turns] = line.split("\t")
                top_cls_turns.append((_relabel_turns(turns), float(energy)))
    return top_cls_turns
