from setuptools import setup, find_packages


def load_requirements(filename):
    with open(filename, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


setup(
    name="quantum-protein-folding-fcc",
    version="0.1.0",
    description="A project for quantum protein folding using FCC lattice models.",
    author="The Cleveland Clinic and IBM Team",
    author_email="your.email@example.com",
    packages=find_packages(),
    install_requires=load_requirements("requirements.txt"),
    python_requires=">=3.10",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
