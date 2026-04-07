from setuptools import setup, find_packages


def load_requirements(filename):
    with open(filename, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


setup(
    name="quantum-protein-folding-fcc",
    version="0.1.0",
    description="A project for quantum protein folding using FCC lattice models.",
    author="The Cleveland Clinic and IBM Team",
    author_email="lir9@ccf.org, raubenb@ccf.org, hakandoga@ibm.com, saki@ibm.com, difilif@ccf.org, radivot@gmail.com, blanked2@ccf.org",
    packages=find_packages(),
    install_requires=load_requirements("requirements.txt"),
    python_requires=">=3.10",
    license="Apache-2.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
