from pathlib import Path
import re

from setuptools import find_packages, setup


REPO_ROOT = Path(__file__).resolve().parent
DATA_INSTALL_PREFIX = "share/quantum-protein-folding-fcc"


def load_requirements(filename):
    with (REPO_ROOT / filename).open(encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def requirement_name(requirement: str) -> str:
    """Return the distribution name for a requirement specifier string."""
    return re.split(r"[<>=!~]", requirement, maxsplit=1)[0].strip()


def collect_data_files(
    source_dir: str, install_prefix: str
) -> list[tuple[str, list[str]]]:
    """Collect non-package data files while preserving directory structure."""
    grouped_files: dict[str, list[str]] = {}
    base_dir = REPO_ROOT / source_dir
    for file_path in base_dir.rglob("*"):
        if file_path.is_file() and file_path.name != ".DS_Store":
            relative_parent = file_path.parent.relative_to(REPO_ROOT)
            target_dir = str(Path(install_prefix) / relative_parent)
            relative_file_path = file_path.relative_to(REPO_ROOT).as_posix()
            grouped_files.setdefault(target_dir, []).append(relative_file_path)
    return [(target_dir, files) for target_dir, files in sorted(grouped_files.items())]


all_requirements = load_requirements("requirements.txt")
optional_requirements = {
    "ray": [requirement for requirement in all_requirements if requirement_name(requirement) == "ray"],
    "runtime": [
        requirement
        for requirement in all_requirements
        if requirement_name(requirement) == "qiskit_ibm_runtime"
    ],
}
install_requires = [
    requirement
    for requirement in all_requirements
    if requirement_name(requirement) not in {"ray", "qiskit_ibm_runtime"}
]
optional_requirements["all"] = optional_requirements["ray"] + optional_requirements["runtime"]


setup(
    name="quantum-protein-folding-fcc",
    version="0.1.0",
    description="A project for quantum protein folding using FCC lattice models.",
    author="The Cleveland Clinic and IBM Team",
    author_email="lir9@ccf.org, raubenb@ccf.org, hakandoga@ibm.com, saki@ibm.com, difilif@ccf.org, radivot@gmail.com, blanked2@ccf.org",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "fcc": ["potentials/*.txt"],
    },
    data_files=[
        (DATA_INSTALL_PREFIX, ["workflow_demo.ipynb"]),
        *collect_data_files("classical_search", DATA_INSTALL_PREFIX),
    ],
    install_requires=install_requires,
    extras_require=optional_requirements,
    python_requires=">=3.10",
    license="Apache-2.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
