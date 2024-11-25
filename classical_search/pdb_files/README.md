Workflow:

1. Download the PDB file of a protein of interest from the Protein Data Bank (PDB) website, e.g.,
    ```
    wget https://files.rcsb.org/download/6VHB.pdb -O "./6VHB.pdb"
    ```
2. Run the `alpha_c.py` script to extract the alpha carbon coordinates from the PDB file and save them in an `xyz` file.
3. Run the classical exhaustive search algorithm (developed by Frank) to obtain the xyz coordinates of the top solutions for FCC and tetrahedral lattices.
4. Use `rmsd` package (`pip install rmsd`) to calculate the RMSD between the experimental structure (alpha carbons only) and the predicted structures, e.g.,
    ```
    calculate_rmsd 6VHB_fcc_1.xyz 6VHB_ca.xyz
    ```
    Note that for this task in we need to set `return_atoms_as_int=False` in lines 1889 and 1896 in the source code `calculate_rmsd.py` to avoid an error.