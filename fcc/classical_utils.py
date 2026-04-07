# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

"""Utility functions relevant for classical search methods."""

import os


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
    new_turns = "".join([mapping[turn] for turn in turns])
    return new_turns[::-1]  # Reverse the string to read from left to right


def _construct_results_file_path(file_name: str, results_dir: str | None = None) -> str:
    """Constructs an absolute path to a classical search results file."""
    if results_dir is None:
        repo_root = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))
        results_dir = os.path.join(repo_root, "classical_search", "results")
    return os.path.normpath(os.path.join(results_dir, file_name))


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
    top_cls_turns = []
    file_path = _construct_results_file_path(file_name, results_dir)
    with open(file_path, "r") as file:
        data = file.read()
        # Split by line
        data_by_line = data.split("\n")
        for line in data_by_line:
            if len(line.split("\t")) == 2:
                [energy, turns] = line.split("\t")
                top_cls_turns.append([_relabel_turns(turns), float(energy)])
    return top_cls_turns
