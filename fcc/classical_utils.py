"""Utility functions relevant for classical search methods."""

import sys

sys.path.append("../")


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


def load_top_cls_solns(file_name: str) -> list[tuple[str, float]]:
    """
    Loads the top solutions from the classical search results file
    (`topojb_<sequence>.txt`).

    Args:
        file_name (str): The name of the file containing the classical search results.

    Returns:
        list[tuple[str, float]]: A list of tuples where each tuple contains a
        string representing the turns and a float representing the energy of
        that configuration.
    """
    top_cls_turns = []
    with open("classical_search/results/" + file_name, "r") as file:
        data = file.read()
        # Split by line
        data_by_line = data.split("\n")
        for line in data_by_line:
            if len(line.split("\t")) == 2:
                [energy, turns] = line.split("\t")
                top_cls_turns.append([_relabel_turns(turns), float(energy)])
    return top_cls_turns
